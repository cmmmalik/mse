from copy import deepcopy
from ase.calculators.vasp import Vasp as asevaspcalc
from ase.atoms import Atoms as aseatoms
from ase.db import connect as asedbconnect
from ase.db.core import Database as dbCore
from ase.db.row import AtomsRow
from ase.io import read as aseread
import numpy as np
import os
import pickle
import warnings

from Database.ase_db import Savedb
from Database.ASE_io import Gpaw_read
from mse.system.directory import CD
from mse.Jobs.job import ASEjob
from Database.ASE_io import Gpawjobdb_read
from mse.analysis.energies import energy_per_atom as simenergy_per_atom, \
    energy_per_formula as simenergy_per_formula
from mse.calc_utilities import get_kpoints_from_density
from mse.system.directory import directorychange
from mse.wrapper.vasp_wrap import generate_runcmd
from mse.utilities import sort_formula_string
from HPCtools.hpc_tools3 import filterargs

# TODO: ions, cell, full relaxation in a single job/job factory
# TODO: create job factory that can do ions, cell, full relaxations


class PickleReadError(Exception):
    pass


class VASP(ASEjob):

    default_files = {"input": "inp.p",
                     "output": "out.p",
                     }

    asedbname = "out.db"

    def __init__(self,
                 name,
                 working_directory: str = None,
                 atoms: aseatoms = None,
                 calcinps: dict = None,
                 modeinps: dict = None,
                 run_type: "static" or "relax" = "static",
                 restart:bool=False):

        super(VASP, self).__init__(name=name,
                                   working_directory=working_directory,
                                   atoms=atoms)
        self._hpc = None
        self._scheme = None
        self._scheme_args = None

        self.run_type = run_type  # depending upon the run_type, fit in appropriate parameters
        self.calcname = "VASP"
        self.inputs = self._default_inputs()
        self.restart = restart
        if restart == True:
            # empty the default outputs...
            self.inputs = {}

        if calcinps is None:
            calcinps = dict()

        self.inputs.update(calcinps)

        if isinstance(modeinps, dict):
            warnings.warn("The use of Mode inps is depracated")
            self.inputs.update(modeinps)

        # outputs
        self._converged = None
        self._row = None
        self._calc = None

    def __repr__(self):
        return super(VASP, self).__repr__()

    def __str__(self):
        return self.__repr__()

    def _default_inputs(self):  # Add default parameters here
        out = dict()
        out["prec"] = "accurate"
        out["ediff"] = 1e-06

        if self.run_type == "relax":
            out.update({"ibrion": 1,
                        "nsw": 60,
                        "ediffg": -1e-02, })

        elif self.run_type == "static":
            out.update({"ibrion": -1,
                        "nsw": 0})

        return out.copy()

    def _directory_change(self):
        cd = CD()
        fpath = os.path.abspath(self.working_directory)
        if fpath != os.getcwd():
            cd.enter(fpath)
        return cd

    @property
    def calc(self):
        return self._calc

    @property
    def vaspcalc(self):
        return self._calc

    @property
    def row(self):
        return self._row
    @row.setter
    def row(self, value):
        if not isinstance(value, AtomsRow):
            raise ValueError("Expected an object of type {}, instead received".format(AtomsRow))
        self._row = value

    @property
    def converged(self):
        return self._converged

    @property
    def scheme(self):
        return self._scheme

    @scheme.setter
    def scheme(self, func):
        if hasattr(func, "__call__"):
            self._scheme = func
        else:
            raise ValueError("Expected a callable function")

    @property
    def scheme_args(self):
        return self._scheme_args

    @scheme_args.setter
    def scheme_args(self, value: dict):
        if isinstance(value, dict):
            self._scheme_args = value
        else:
            raise ValueError("must be an instance of {}".format(type(value)))

    @property
    def run_type(self):
        return self._run_type

    @run_type.setter
    def run_type(self, value: "static" or "relax"):
        if value in ["static", "relax"]:
            self._run_type = value
        else:
            raise ValueError("'run_type' must be 'static' or 'relax', got value as '{}'".format(value))

    def run(self):
        if not self.scheme:
            out = self.run_conventional()
        else:
            out = self.run_scheme()
        self._setfinished()
        return out

    def run_conventional(self):
        return self.atoms.get_potential_energy()

    def run_scheme(self):
        torun = self.intialize_scheme()
        if hasattr(torun, "run"):
            return torun.run()
        return torun

    def intialize_scheme(self):
        #TODO: need thorough testing
        torun = self.scheme(atoms=self.atoms, **self.scheme_args)
        return torun

    def initialize_calc(self, **kwargs):

        if self.atoms and getattr(self.atoms,"calc"):
            return

        self.inputs.update(kwargs)
        self.set_kden() # check for kpts

        calc_inputs = deepcopy(self.inputs)
        calc_inputs.pop("kpden", None)

        if not self.restart: # let vasp calculator handle the job reading...

            calc = asevaspcalc(atoms=self.atoms,
                               **calc_inputs)

        else:
            calc, atoms = self.readcalc(**calc_inputs)
            self.atoms = atoms

        self.atoms.calc = calc

    def readcalc(self, **calc_inputs):

        calc = asevaspcalc(restart=True, **calc_inputs)
        atoms = calc.get_atoms()
        calc.reset()
        return calc, atoms

    def initialize(self):
        if not os.path.exists(self.working_directory):
            super(ASEjob, self).initialize()

        cd = self._directory_change()
        self._save()
        cd.exit()
        self._setinitialized()

    def _save(self):
        #TODO: find a better way to pickle the object, may be with ase.ulm
        with open(self.default_files["input"], "wb") as f:
            pickle.dump(self, f)

    def set_kden(self):
        kden = self.inputs.get("kpts")
        if isinstance(kden, dict):
            if "density" in kden:
                kpts = get_kpoints_from_density(kden["density"],
                                                cell=self.atoms.get_cell(complete=True),
                                                even=kden.get("even", False))
                self.inputs.update({"kpts": kpts,
                                    "kpden": kden})  # do not loose the original input


    @classmethod
    @directorychange
    def read_job(cls, filename=None, verbosity: int = 0):
        if not filename:
            filename = cls.default_files["output"]
        try:
            with open(filename, "rb") as f:
                job = pickle.load(f)
        except Exception as error:
    #       warnings.warn("Could not read the {}, either it does not exist or un-readable, "
    #                      "reading input file {}".format(filename, cls.default_files["input"]))
            raise PickleReadError("Unable to unpickle the job") from error

        # read_row
        try:
            row = cls.ase_row_db(db=cls.asedbname)
            atoms = row.toatoms(False)
            print("Atoms read from the ase row")
        except FileNotFoundError:
            warnings.warn("Ase db file '{}' does not exist".format(cls.asedbname))
            warnings.warn("Either ask calculator for the output parameters i.e. energies, forces or, if present,"
                  ".txt file will be read ")
            row = None
            try:
                atoms = aseread("OUTPUT")  # reading only the last configuration
                print("Atoms read from OUTPUT file")
            except (IOError, OSError, IndexError) as ex:
                warnings.warn("Invalid/empty file: '{}'".format("OUTPUT"))
                raise ex

        job.atoms = atoms
        job.row = row

 #       if not job.vaspcalc:
 #           calc = asevaspcalc(restart=True)
 #           job._calc = calc

        return job

    @classmethod
    @directorychange
    def read_row(cls, verbosity: int = 0):
        try:
            row = cls.ase_row_db(cls.asedbname, verbosity=verbosity)
        except FileNotFoundError:
            warnings.warn("Ase db file '{}' does not exist".format(cls.asedbname), RuntimeWarning)
            row = None
        return row

    @classmethod
    @directorychange
    def read_atoms(cls):
        try:
            atoms = aseread("OUTCAR")
        except OSError:
            warnings.warn("Invalid/empty file: OUTCAR")
            atoms = None
        return atoms

    def _savetodb(self):
        with asedbconnect(self.asedbname) as mydb:
            mydb.write(self.atoms)

    def add_to_ase_database(self, db: str or dbCore, **external_keys):
        keys, data = self._vaspdb_keys()
        row = self.row

        if not row:
            row = self.atoms

        formula = self.atoms.get_chemical_formula(empirical=True)
        if keys["path"]:
            keys["rpath"] = os.path.relpath(keys["path"])

        keys["wrkdir"] = "{}".format(self.working_directory)
        keys["static"] = self.run_type == "static"
        keys["relaxed"] = not keys["static"]
        keys["primitive_formula"] = sort_formula_string(formula, )

        keys.update(external_keys)

        print(row)
        self._smartwriteasedb(db, row=row, data=data, keys=keys)

    @directorychange
    def _vaspdb_keys(self):
        keys = dict()
        data = dict()

        row = self.row
        subdir = self.inputs.get("directory", ".")

        try:
            params = row.calculator_parameters.copy()
            asepath = os.path.abspath(self.asedbname)

        except AttributeError:
            try:
                params = self._vasp_calc_parameters()
            except AttributeError:
                if self.status != "finalized":
                    raise RuntimeError("Job is not finalized, cannot used 'inputs'")
                params = self.inputs.copy()


            asepath = None
        encut = params.get("encut")
        if asepath:
            keys["asepath"] = asepath
        jobpath = os.path.abspath(self.default_files["output"])
        if os.path.exists(jobpath):
            keys["jobpath"] = jobpath
        path = os.path.abspath(os.path.join(subdir, "OUTCAR"))
        if os.path.exists(path):
            keys["path"] = path
        if encut:
            keys["encut"] = encut

        keys.update(Gpawjobdb_read.get_keys_asedb(params))
        kkpts = self._get_kden_kpts()

        if kkpts.get("kpden"):
            gamma = params.get("gamma")
            if gamma:
                kkpts["kpden"]["gamma"] = params.get("gamma")

        keys.update(kkpts)
        data.update(kkpts)

        if "kpts" not in keys:  # ask the calculator
            keys["kpts"] = self.calc.kpts

        # energy_per_formula and energy_per_atom
        keys.update({"energy_per_atom": simenergy_per_atom(row=row),
                     "energy_per_formula": simenergy_per_formula(row=row)})

        init_poscar = Savedb.find_poscar(directory=subdir, query=["POSCAR"])
        if isinstance(init_poscar, (list, tuple)):
            init_poscar.sort()
            init_poscar = init_poscar[0]

        if not init_poscar:  # highly unlikely
            init_poscar = os.path.join(subdir, "OUTCAR")
        try:
            init_atoms = aseread(init_poscar, index=0, format="vasp")
        except OSError:
            warnings.warn("Couldnt' read the initial poscar")
            init_atoms = None

        if init_atoms:
            keys["init_poscar"] = True
            keys["poscarname"] = os.path.abspath(init_poscar)

            constraints = init_atoms.constraints #check for constraints
            if constraints and isinstance(constraints, (list, tuple)):
                    constraints = [i.todict() for i in constraints]
                    data["init_atoms_constraints"] = constraints # save into dictionary
                    del init_atoms.constraints # remove constraints

            data["init_atoms"] = init_atoms
        keys = Gpaw_read.keys_parsedb(inkeys=keys)
        return keys, data

    def _get_kden_kpts(self):
        keys = {}
        kden = self.inputs.get("kpden", [])
        kpts = self.inputs.get("kpts", [])
        if np.asarray(kden).size:
            keys["kpden"] = kden
            keys["kden"] = kden["density"]
        if np.asarray(kpts).size:
            keys["kpts"] = kpts

        return keys

    def _vasp_calc_parameters(self):
        params = {}
        calc = self.calc
        if not calc:
            raise RuntimeError("Calc is not present, cannot read the parameters")

        for kk, vv in calc.param_state.items():
            for k, v in vv.items():
                if v is not None:
                    if isinstance(v, np.ndarray):
                        if v.size == 0:
                            continue
                    elif v is not False and not v:
                        continue
                    params[k] = v
        params.update({"txt": calc.txt})
        return params

    def submit_hpc(self, dry_run: bool = True,
                   save_calc: bool = True,
                   save_db: bool = True,
                   backup: bool = False,
                   remove: bool = True,
                   envcmd: bool = None,
                   module: str = None,
                   code: str = None,
                   **kwargs):

        self.run_check()
        pre_runcmd = kwargs.pop("pre_runcmd", None)
        submitargs, prepargs = filterargs(self.hpc.server, **kwargs)
        print("submission arguments : {}".format(submitargs))
        print("preparation arguments: {}".format(prepargs))

        if not "runcmd" in submitargs:

            submitargs["runcmd"] = generate_runcmd(inputfile=self.default_files["input"],
                                                   save_calc=save_calc,
                                                   save_db=save_db,
                                                   backup=backup,
                                                   remove=remove,
                                                   envcmd=envcmd)
        if not module and not code:
            code = "vasp"


        if pre_runcmd:
            submitargs["runcmd"] = pre_runcmd + "\n" + submitargs["runcmd"]

        self.hpc.server.write_submit(code=code, module=module, **submitargs)

        if not dry_run:
            self.hpc.prepare(**prepargs)
            self.hpc.submit()

    def todict(self):
        dct = {}
        for k, v in self.__dict__.items():
            if isinstance(v, aseatoms):
                v = v.todict()

            elif k == "_scheme":
                continue

            elif hasattr(v, "todict"):
                v = v.todict()

            elif hasattr(v, "asdict"):
                warnings.warn("asdict() method of {} if undepracated",DeprecationWarning)
                v = v.asdict()

            elif k == "_working_directory":
                v = str(v)

            else:
                print("Can not write {} instance of {} as {}".format(k,v,dict))

            dct[k] = v

        return dct

    def as_dict(self):
        dct = {}
        dct["@module"] = self.__class__.__module__
        dct["@class"] = self.__class__.__name__
        dct.update(self.todict())
        return dct

    @classmethod
    def from_dict(cls, dct):
        from monty.json import MontyDecoder
        dct["inputs"] = MontyDecoder().process_decoded(dct["inputs"])

        obj = cls(name=dct["_name"],
                  working_directory=dct["_working_directory"],
                  calcinps=dct["inputs"])

        for k,v in dct.items():
            if k in ["_name", "_working_directory"] or k.startswith("@"):
                continue
            elif k == "_atoms":
                v = aseatoms.fromdict(MontyDecoder().process_decoded(v))
            setattr(obj, k, v)

        return obj

    def tojson(self):

        from monty.json import MontyEncoder
        import json
        return json.dumps(self, cls=MontyEncoder)

    @classmethod
    def fromjson(cls, jstring):
        from monty.json import MontyDecoder
        import json

        dct = json.loads(jstring, cls=MontyDecoder)
        return cls.from_dict(dct)
