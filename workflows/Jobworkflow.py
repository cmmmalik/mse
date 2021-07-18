import os
import sys

from ase import Atoms
from ase.io import read as aseread
from ase.db.core import Database as DBcore
from ase.db import connect
from ase.parallel import paropen, parprint
import numpy as np
import json
import pickle

from mse.calculator.gpaw_calc import Gpaw
from mse.wrapper.gpaw_wrap import initialize_calc, savetodb
from mse.system.directory import directorychange
from mse.analysis.energies import energy_per_atom as get_energy_per_atom, energy_per_formula as get_energy_per_formula
from Database.ASE_io import Gpawjobdb_read
from Database.ase_db import Savedb


class Smearingworkflow:

    """
    A basic workflow that will perform smearing study using a calculator job i.e. gpaw or vasp (at the moment).
    The main usage is for doing the convergence w.r.t smearing. This is a very basic implementation.
    ToDo: Need to generalise this and make it efficient.

    Basic Usage and Examples:
    """
    def __init__(self,
                 name: str,
                 working_directory: str,
                 atoms: Atoms,
                 calc_type: "gpaw" or "vasp" = "gpaw",
                 smearing_name: str ="fermi",
                 widthlst: list or tuple = None,
                 widthlimits: [int, int] = None,
                 calcinps: dict = None,
                 modeinps: dict = None,
                 verbosity: int=1,
                 **kwargs,
                 ):

        assert widthlst or widthlimits
        self._job = None
        self._smearing = None
        self.smearingname = smearing_name
        self._logfile = "output.log"
        self._atoms = None
        self._outputs = dict()

        if isinstance(widthlst, (list, tuple)):
            self._widths = widthlst
        else:
            if not len(widthlimits) == 2:
                raise ValueError("Expected a list or tuple of two values, but got {}".format(widthlimits))
            self._widths = np.arange(widthlimits[0], widthlimits[-1])

        if calc_type != "gpaw":
            raise NotImplementedError("Other than 'gpaw', not implemented yet")
     #       job = Gpaw(name=name,
     #                  working_directory=working_directory,
     #                  atoms=atoms,
     #                  run_type="static",
     #                  calcinps=calcinps,
     #                  modeinps=modeinps,
     #                  **kwargs)

        self._calc_type = calc_type

        self.verbosity = verbosity
        self.atoms = atoms

        # job attributes, that will be used to initialze and run the job
        jobattributes = {}
        jobattributes["name"] = name
        jobattributes["working_directory"] = working_directory
        jobattributes["calcinps"] = calcinps
        jobattributes["modeinps"] = modeinps
        jobattributes.update(kwargs)

        jobattributes["run_type"]="static"

        self._jattr = jobattributes

    @property
    def smearing_name(self):
        return self._smearing_name

    @smearing_name.setter
    def smearing_name(self, name:str):
        if name in ["fermi"]:
            self._smearing = name

    def todict(self):
        dct = {}
        for k,v in self.__dict__.items():
            if isinstance(v, (Gpaw)): # I am loosing the job here
                continue
            elif isinstance(v, Atoms):
                v = v.todict()
            dct[k] = v

        return dct

    def tojson(self, *args, **kwargs):
        return json.dumps(self.todict(), *args, **kwargs)

    @property
    def job(self):
        return self._job

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, value:Atoms):
        if not isinstance(value, Atoms):
            raise ValueError("Expected an {} object instead got {}".format(Atoms.__name__,type(value)))
        self._atoms = value

    @property
    def energies(self):
        return self._outputs.get("energies", None)

    def _setup_smearing(self):
        if self.verbosity >= 1:
            print("Removing occupations from the inputs, if already present")

        if isinstance(self.job, Gpaw):
            self.job.inputs.pop("occupations", None)

    def set_job(self):
        assert self.atoms
        if self._calc_type == "gpaw":
            job = Gpaw(**self._jattr)
            self._job = job

    def setup(self):
        self.set_job()
        self._setup_smearing()
        initialize_calc(job=self._job) # will initialize the actual GPAW calculator with PW mode
        sst = " {:>12} {:>12} {:>12}".format("No.", "width", "energy")
        with paropen(self._logfile, "a+") as ff:
            ff.write(sst)

    def run(self):
        job = self.job
        calc = job.atoms.calc

        for i, width in enumerate(self._widths):
            if self.verbosity >=1:
                parprint('Iteration No. : {}'.format(i), flush=True)  # only master
            calc.set(occupations={"name": self.smearingname, "width": width})
            job.static_run() # ask calculator of the job to do the calculation.
            ### here do the logging, before moving forward
            with paropen("output.log", "a+") as ff:
                ff.write("{:>12} {:>12} {:>12}".format(i, width, job.atoms.get_potential_energy()))
            #add parameters
            keys = dict()
            keys["smearing_name"] = self.smearingname
            keys["width"] = width
            keys.update(self._genkeys(params=calc.parameters.copy()))
            data, _keys = self._gendata()
            keys.update(_keys)
            keys.update({"energy_per_atom" : get_energy_per_atom(atoms=job.atoms),
                        "energy_per_formula": get_energy_per_formula(atoms=job.atoms)})

            savetodb(job, data=data, **keys) #save to ase db row on master node only


    def _genkeys(self, params: dict ):
        keys = {}
        mode = params.get("mode")
        keys["encut"] = mode.get("ecut") if isinstance(mode, dict) else mode.todict().get("ecut")
        keys.update(Gpawjobdb_read.get_keys_asedb(params))

        try:
            keys.update(self.job._addkpdens_kpts(keys=keys))
        except AttributeError:
            pass
        return keys


    def _gendata(self):

        data = dict()
        keys = dict()

        init_poscar = Savedb.find_poscar(directory=".")
        if isinstance(init_poscar, list):
            init_poscar = init_poscar[0]
        if not init_poscar:
            init_poscar = self.job.inputs["calc_args"]["txt"]
        try:
            init_atoms = aseread(init_poscar, index=0)
        except OSError:
            init_atoms = None
        #fill in
        if init_atoms:
            data["init_atoms"] = True
            keys["poscarname"] = init_poscar
            keys["init_poscar"] = True

        return data, keys


    def _finalize(self):
        self.job = None
        with paropen("out.p", "wb") as f:
            pickle.dump(self, f)
        with open("out.json") as ff:
            json.dump(self.todict(), ff)

    @staticmethod
    def fromdict(dct):
        copydct = dct.copy()
        atoms = copydct.pop("_atoms")
        atoms = Atoms.fromdict(atoms)
        _outputs = copydct.pop("_outputs", None)
        _jattrs = copydct.pop("_jattr", None)
        run_type = _jattrs.pop("run_type", None)

        obj = Smearingworkflow(atoms=atoms, **copydct, **_jattrs)
        if run_type:
            obj._jattr["run_type"] = run_type
        if _outputs:
            obj._outputs = _outputs
        return obj

    @classmethod
    @directorychange
    def read(cls, filename=None, use_pickle=True):
        if not filename:
            filename = "out.p"
        unable_topick = False

        if use_pickle:
            try:
                with open(filename, "rb") as f:
                    obj = pickle.load(f)
            except OSError:
                print("Can't pickle the job, looking for .json file")
                unable_topick = True

        if not use_pickle or unable_topick:
            with open("{}.json".format(filename.rsplit(".", maxsplit=1)[0]), "rb") as jf:
                dct = pickle.load(jf)
                obj = cls.fromdict(dct)

        return obj

    def read_rows(self):
        paths = os.path.join(self._jattr["working_directory"],"out.db")
        with connect(paths) as db:
            rows = [r for r in db.select()]
        return rows

    def add_to_ase_database(self, db: str or DBcore,add_extra_keys: bool = True, **externalkeys):
        rows = self.read_rows()

        if not rows:
            raise ValueError("Unable to read rows, can't proceed further")
        if add_extra_keys:
            keys = self._rowget_keys()
        else:
            keys = dict()

        if isinstance(db, str):
            db = connect(db)
        try:
            db.__enter__()
            for r in rows:
                db.write(r, key_value_pairs=keys, **externalkeys)
        finally:
            exc_type, exc_value, tb = sys.exc_info()
            db.__exit__(exc_type, exc_value, tb)

# safe exit

    def _rowget_keys(self):
        keys = dict()
        #paths here
        asepath = os.path.join(self._jattr["working_directory"], "out.db")
        jobpath = os.path.join(self._jattr["working_directory"], "out.p")

        if os.path.exists(jobpath):
            keys["jobpath"] = jobpath

        if os.path.exists(asepath):
            keys["asepath"] = asepath
        keys["working_directory"] = self._jattr["working_directory"]
        return keys

