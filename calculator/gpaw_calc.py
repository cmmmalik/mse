from ase.db import connect as asedbconnect
from ase.db.core import Database as dbCore
from ase.io import read as aseread
from ase.io.trajectory import Trajectory
import pickle
from gpaw import GPAW as gpawGPAW
import os
import warnings

from Database.ase_db import Savedb
from Database.ASE_io import Gpawjobdb_read, Gpaw_read
from mse.Jobs.Gpaw import Gpaw as Gpawjob
from mse.wrapper.gpaw_wrap import generate_runcmd
from mse.optimizer.optimizer import Optimize, Moptimize
from mse.system.directory import CD
from mse.io.gpaw_io import Readparameters
from mse.analysis.energies import energy_per_atom as simenergy_per_atom, \
    energy_per_formula as simenergy_per_formula
from HPCtools_v2.hpc_tools3 import filterargs, directorychange

# TODO: add relaxation schemes
# TODO: Add save_to_database function with columns keys of parameters. (See Gpawjob class and ase_db/ASE_IO)


class Gpaw(Gpawjob):

    # defaults_inputs = {"calc_args": {"kpts": {"density": 4, "gamma": True},
    #                    "xc": "PBE",
    #                    "eigensolver": "rmm-diis",
    #                    "occupations": {"name": "fermi-dirac", "width": 0.4}},
    #                    "mode_args": {"encut": 500}}
    # This is not a good approach for keeping a dictionary here like this
    # Use mutable objects like tuple, etc.

    defaults_files = {"input": "inp.p",
                      "output": "out.p",
                      "calc": "calc.gpw",
                      "traj": "out.traj"} # The things whill will not be changes
    # TODO: Still change them into a tuple

    #outputs
    asedbname = "out.db"

    def __init__(self,
                 name,
                 working_directory=None,
                 atoms=None,
                 run_type: "static" or "relax" = "static",
                 calcinps: dict = None,
                 modeinps: dict = None,
                 relaxinps: dict = None,
                 **kwargs):

        super(Gpaw, self).__init__(name=name, working_directory=working_directory, atoms=atoms)
        self.calc = "GPAW"
        self.mode = "PW"
        self.inputs = self._default_inputs()  # Always copy the dictionary
        self.run_type = run_type
        self._optimizer = None
        self._newrunscheme = None

        if self.run_type == "static":
            self.inputs["calc_args"]["txt"] = "static.txt"
            if relaxinps:
                warnings.warn("{} was provided but will not be used, as run_type {} was chosen".format(relaxinps,
                                                                                                       run_type))
        elif self.run_type == "relax":
            if relaxinps is None:
                relaxinps = dict()
            self._initrelax_inputs(**relaxinps)

            # set trajectory
            if kwargs.get("trajectory", False):
                self.relax_inputs["trajectory"] = self.defaults_files["traj"]
        else:
            raise ValueError("run_type can have 'static' or 'relax' values, but got {}".format(self.run_type))

        if isinstance(calcinps, dict):
            if calcinps is None:
                calcinps = dict()

            self.inputs["calc_args"].update(calcinps)

        elif calcinps is not None:
            raise TypeError("'calcinps' must be an instance of {}, but found an instance of type {} ".format(dict,
                                                                                                             type(calcinps)))

        if isinstance(modeinps, dict):
            if modeinps is None:
                modeinps = dict()
            self.inputs["mode_args"].update(modeinps)

        elif modeinps is not None:
            raise TypeError("'modeinps' must be an instance of {}, but found an instance of type {}".format(dict,
                                                                                                            type(modeinps)))

        # outputs
        self.converged = None
        self.traj = None
        self.row = None
        self.scf = None

    def _initrelax_inputs(self, **kwargs):

        self.inputs["calc_args"]["txt"] = "rel.txt"
        self.relax_inputs = self._relax_defaults()
        self.relax_inputs.update(kwargs)

    @staticmethod
    def _default_inputs():

        out = {"calc_args": {"kpts": {"density": 4, "gamma": True},
                                         "xc": "PBE",
                                         "eigensolver": "rmm-diis",
                                         "occupations": {"name": "fermi-dirac", "width": 0.4}},
                           "mode_args": {"encut": 500}}
        return out

    @staticmethod
    def _relax_defaults():
        return {"fmax": 0.01,
                "reltype": "ions" }

    def static_run(self):
        self.atoms.get_potential_energy()

    def relax_run(self):
        self.setoptimizer()
        return self.optimizer.optimize()

  #  def relax_cell(self): # run should call this
 #       optimize(atoms=self.atoms, reltype="cell", fmax=self.fmax)

  #  def relax_all(self):
  #      optimize()

    def run(self):

        if self.run_type == "static":
            self.static_run()
        elif self.run_type == "relax":
            return self.relax_run()
        else:
            raise ValueError("Unknown value of 'run_type': {}".format(self.run_type))

        self._setfinished()

    def submit_hpc(self, dry_run: bool = True,
                   save_calc: bool = True,
                   save_db: bool = False,
                   backup: bool = False,
                   remove: bool = True,
                   envcmd: bool = None,
                   module: str = None,
                   code: str = None,
                   **kwargs):

        self.run_check()
        submitargs, prepargs = filterargs(self.hpc.server, **kwargs)
        if not "runcmd" in submitargs:

            submitargs["runcmd"] = generate_runcmd(inputfile=self.defaults_files["input"],
                                                   save_calc=save_calc,
                                                   save_db=save_db,
                                                   backup=backup,
                                                   remove=remove,
                                                   envcmd=envcmd)
        if not module and not code:
            code = "gpaw"

        self.hpc.server.write_submit(code=code, module=module, **submitargs)

        if not dry_run:
            self.hpc.prepare(**prepargs)
            self.hpc.submit()

    @classmethod
    @directorychange
    def read_job(cls, filename=None, attach_calc: bool = False, trajectory: str = None):
        atoms = None

        if not filename:
            filename = cls.defaults_files["output"]

        with open(filename, "rb") as f:
            job = pickle.load(f)

        try:
            filename = cls.defaults_files["calc"]
            calc = gpawGPAW(filename)
        except FileNotFoundError:
            warnings.warn(f"File {filename} does not exist")
            calc = None
        except Exception as e:
            print("**Encountered Unexpected exception {}".format(e))
            calc = None

        # read db rows
        try:
            row = cls.ase_row_db()
            atoms = row.toatoms(attach_calc) # TODO: need to do this in a better way
            print("Atoms read from the ase row")
            calc = atoms.calc
            # TODO: Currently, ase db stores only the last configuration either shift
            #       to trajectory (done) or stores all the configurations in the database
            #       (advantage will be of saving the calculator along with initial and final parameters).
            #       A good alternative to saving the big  gpw-format file of the calculator

        except FileNotFoundError:
            print("Ase db file '{}' does not exist".format(cls.asedbname))
            print("Either ask calculator for the output parameters i.e. energies, forces or, if present,"
                  ".txt file will be read ")
            row = None
            try:
                atoms = aseread(job["calc_args"]["output"])  # reading only the last configuration
                print("Atoms read from .txt file")
            except IOError:
                print("Invalid/empty file: '{}'".format(job["calc_args"]["output"]))
                pass

        #read trajectory_file
        try:
            traj = cls.read_traj(filename=trajectory)
            job.traj = traj

        except IOError:
            warnings.warn("Couldn't read the trajectory file {}".format(cls.defaults_files["traj"]))

        if atoms:
            job.atoms = atoms
        else:
            if job.atoms:
                warnings.warn("Job(unpickled) already contains the atoms"
                              "since atoms object couldn't read from the directory. "
                              "warning can be ignored. Just to be on safe side, inspect the atoms")
        job.row = row

        return job, calc

    @classmethod
    def read_traj(cls, filename: str = None):

        if not filename:
            filename = cls.defaults_files["traj"]

        return Trajectory(filename=filename)

    @classmethod
    def read_db(cls):
        atoms = aseread(cls.asedbname)
        return atoms

    @classmethod
    def ase_row_db(cls, index: int = -1):
        rows = Gpawjobdb_read(db=cls.asedbname).R
        if len(rows) > 1:
            warnings.warn("Found more than one rows in the database {}"
                          "\nSelecting the last index row\nUse 'index' to change the read row".format(cls.asedbname))
        rows = rows[index]
        return rows

    @classmethod
    def read_atoms_db(cls):
        return cls.ase_row_db().toatoms(True)

    def _save(self):
        with open(self.defaults_files["input"], "wb") as f:
            pickle.dump(self, f)

    def _savetodb(self):
        with asedbconnect(self.asedbname) as mydb:
            mydb.write(self.atoms)
        print("Successfully written into the database")

    def _dirchanger(self):
        cd = CD()
        if self.working_directory != os.getcwd():
            cd.enter(self.working_directory)
        return cd

    def initialize(self):
        if os.path.exists(self.working_directory):
            super(Gpawjob, self).structure_write(filename=os.path.join(self.working_directory,
                                                getattr(self.inputs, "poscarname", "POSCAR")),
                                                 format="vasp",
                                                 vasp5=True)
        else:
            super(Gpaw, self).initialize()

        cd = self._dirchanger() # Jump in
        self._save()
        cd.exit() # Jump out
        self._setinitialized()

    def setoptimizer(self):
        if not self.newrunscheme:
            self.optimizer = Optimize(atoms=self.atoms, **self.relax_inputs)
        else:
            self.optimizer = self.newrunscheme(atoms=self.atoms, **self.relax_inputs)

    @property
    def optimizer(self):
        return self._optimizer

    @optimizer.setter
    def optimizer(self, value: Optimize or Moptimize):
        if isinstance(value, Optimize) or isinstance(value, Moptimize):
            self._optimizer = value
        else:
            raise TypeError("Expected an instance of {0}, but got {1}".format(Optimize, type(value)))

    @property
    def newrunscheme(self):
        return self._newrunscheme

    @newrunscheme.setter
    def newrunscheme(self, value):
        self._newrunscheme = value

    @classmethod
    def help(cls):
        """
        help method to display the usage of a class
        :return:
        """
        pass

    def add_to_ase_database(self, db: str or dbCore, **externalkeys):

        outdb = os.path.join(self.working_directory, self.asedbname)

        if self.row or os.path.exists(outdb):
             keys, data = self._newgpaw_keys()
             row = self.row
        else:
            gatoms, data, keys = self._oldgpaw_reader()

            if gatoms != self.atoms:
                warnings.warn(f"Atoms read by Gpaw reader and already-present in the instance are not same", ValueError)
            row = gatoms

        if not "rpath" in keys:
            keys["rpath"] = os.path.relpath(keys["path"])

        keys["wrkdir"] = "{}".format(self.working_directory)
        keys["static"] = self.run_type == "static"
        keys["relaxed"] = not keys["static"]
        keys.update(externalkeys)
        print("Row:\n{}".format(row))
        self._smartwriteasedb(db, row=row, data=data, keys=keys)

    def _smartwriteasedb(self, db: str or dbCore, row, data, keys):
        try:
            db.write(row, data=data, **keys)
        except (AttributeError, IOError) as e:
            with asedbconnect(db) as mydb:
                mydb.write(row, data=data, **keys)

    def _oldgpaw_reader(self):

        gpaw_reader = Savedb.gpawcollect(txtfile=self.inputs["calc_args"]["txt"],
                                         pyfile=None,
                                         poscarfile=None,
                                         directory=self.working_directory)

        atoms = gpaw_reader.atoms
        gpaw_reader.get_keys_asedb(kpts={})

        data = gpaw_reader.output_data_asedb()
        data.update(self.inputs)
        keys = gpaw_reader.keys_parsedb(inkeys=data)

        return atoms, data, keys

    @directorychange
    def _newgpaw_keys(self):

#TODO: function is too big, split it

        keys = dict()
        row = self.row
        if not row:
            try:
                row = self.ase_row_db()
                self.row = row
            except IOError:
                pass

        try:
            params = row.calculator_parameters.copy()
            asepath = os.path.abspath(self.asedbname)
            encut = params["mode"].get("ecut", None)
        except AttributeError:
            if self.status != "finalized":
                raise RuntimeError("Job is not finalized, cannot used 'inputs'")

            params = self.inputs.copy()
            encut = params["mode_args"].get("encut", None)
            asepath = None

        if asepath:
            keys["asepath"] = asepath

        jobpath = os.path.abspath(self.defaults_files["output"])
        if os.path.exists(jobpath):
            keys["jobpath"] = jobpath

        path = os.path.abspath(self.inputs["calc_args"]["txt"])
        if os.path.exists(path):
            keys["path"] = path

        if encut:
            keys["encut"] = encut

        keys.update(Gpawjobdb_read.get_keys_asedb(params))
        keys.update(self._addkpdens_kpts(keys=keys))

        if not "kpts" in keys:
            # read from .txt file
            calc_inparameters = Readparameters(file=self.inputs["calc_args"]["txt"])
            try:
                keys["kpts"] = calc_inparameters.parameters["k-points"]
            except KeyError:
                pass

        init_poscar = Savedb.find_poscar(directory=".")
        if isinstance(init_poscar, list):
            init_poscar = init_poscar[0]

        if not init_poscar:
            init_poscar = self.inputs["calc_args"]["txt"]
        try:
            init_atoms = aseread(init_poscar, index=0)
        except OSError:
            init_atoms = None

        data = dict()
        if init_atoms:
            keys["init_poscar"] = True
            keys["poscarname"] = os.path.abspath(self.inputs["calc_args"]["txt"])
            data["init_atoms"] = init_atoms
        #add energies per formuls

        keys.update({"energy_per_atom": simenergy_per_atom(row=row),
             "energy_per_formula": simenergy_per_formula(row=row)})

        keys = Gpaw_read.keys_parsedb(inkeys=keys)
        return keys, data

    @staticmethod
    def _addkpdens_kpts(keys: dict):

        kpts = keys.pop("kpts")

        out = {}
        if not kpts:
            return out # My work is done here, immediate leave

        if "density" in kpts:
            out["kpden"] = kpts
            out["kden"] = kpts["density"]

        elif "size" in kpts:
            out["kpts"] = kpts

        return out