from ase.atoms import Atoms
from ase.spacegroup.symmetrize import refine_symmetry, check_symmetry
from collections.abc import Iterable
from copy import deepcopy
import numpy as np
import os
import warnings

from mse.calculator.gpaw_calc import Gpaw
from mse.calculator.vasp_calc import VASP
from mse.Jobs.job import HPCjob
from HPCtools.hpc_tools3 import filterargs
# TODO: Implement proper Workflow classes (Next steps); Below is just a proof of concept


class Workflow:

    def __init__(self,
                 encut: int or float,
                 atoms: Atoms = None,
                 kpts: list or tuple or Iterable = None,
                 kplimits: [int, int] = None,
                 calculator_type: str = "gpaw",
                 verbosity: int = 1):

        assert kpts or kplimits

        if isinstance(kpts, Iterable):
            self._kpts = kpts
        else:
            if not len(kplimits) == 2:
                raise ValueError("Expected a list or tuple of two values, but got {}".format(kplimits))
            self._kpts = np.arange(kplimits[0], kplimits[-1])

        self._encut = encut
        self._atoms = atoms
        self._calculator_type = calculator_type
        self.verbosity = verbosity

        self._wrkpaths = None
        self._jobs = None

    @property
    def atoms(self):
        return self._atoms

    @property
    def encut(self):
        return self._encut

    @property
    def kpts(self):
        return self._kpts

    @property
    def calculator_type(self):
        return self._calculator_type

    @property
    def wrkpaths(self):
        return self._wrkpaths
    @property
    def jobs(self):
        return self._jobs

    def generate_working_paths(self, outerdirectory: str = "."):

        encut = self.encut
        if not isinstance(encut, Iterable):
            encut = [encut]

        wrkpaths = [os.path.join(outerdirectory, "{}.{}").format(en, k) for en in encut for k in self.kpts]

        if self.verbosity >= 1:
            print("Workpaths (subdirectories): {}".format(wrkpaths))

        self._wrkpaths = wrkpaths
        return wrkpaths

    def check_symmetry(self, symprec=1e-06):
        return check_symmetry(self.atoms, symprec=symprec)

    def refine_symmetry(self, symprec=0.01):
        refine_symmetry(self.atoms, symprec=symprec)

    def pre_setup_jobs(self, outerdirectory: str = "", symmetrize: bool = True, symprec=1e-03):
        self.generate_working_paths(outerdirectory=outerdirectory)
        if symmetrize:
            self.refine_symmetry(symprec=symprec)

    def get_kpts(self, kpts: dict or tuple = ("gamma", True)):

        if not isinstance(kpts, dict):
            try:
                kpts = dict(kpts)
            except (TypeError, ValueError) as error:
                kpts = dict((kpts,))
        if self.verbosity >= 2:
            print("Debug:\nkpts: {}".format(kpts))
        Kpts = []
        for kden in self.kpts:
            ckpts = deepcopy(kpts)
            ckpts.update({"density": kden})
            Kpts.append(ckpts)

        return Kpts

    def setup_jobs(self,
                   name: str,
                   kpts: dict or tuple = ("gamma", True),
                   modeinps: dict = None,
                   calcinps: dict = None,
                   run_type: 'static' or 'relax' = "static",
                   **kwargs):

        if not calcinps:
            calcinps = dict()
        if not modeinps:
            modeinps = dict()

        Kpts = self.get_kpts(kpts=kpts)
        if self.verbosity >= 2:
            print("Debug:\nKpts: {}".format(Kpts))

        if self.calculator_type == "gpaw":
            Jobfunc = Gpaw
        elif self.calculator_type == "vasp":
            Jobfunc = VASP
            warnings.warn("Mode inps:{}, will not be used".format(modeinps))
            print("For VASP they are part of calcinps")
            print("Cleaning....")
            modeinps = dict()

        else:
            raise ValueError("Unknown {} calculator_type given".format(self.calculator_type))

        jobs = []
        for i, w in enumerate(self.wrkpaths):
            c_calcinps = deepcopy(calcinps)
            c_calcinps.update({"kpts": Kpts[i]})

            c_modeinps = deepcopy(modeinps)
            c_modeinps.update({"encut": self.encut})

            j = Jobfunc(name=name,
                     working_directory=w,
                     atoms=self.atoms,
                     run_type=run_type,
                     modeinps=c_modeinps,
                     calcinps=c_calcinps,
                     **kwargs)

            jobs.append(j)

        self._jobs = jobs

        return jobs

    def _show_attr(self, attrs):
        st = ["{}\n".format(getattr(j, attrs)) for j in self.jobs]
        st = "\n".join(st)
        st = "{}:\n{}".format(attrs, st)
        return st

    def show_inputs(self):
        st = self._show_attr("inputs")
        print(st)

    def _auto_fetch(self, value: bool = True):
        list(map(lambda j: setattr(j, "autofetch", value), self.jobs))

    def _initialize_jobs(self):
        list(map(lambda j: j.initialize(), self.jobs))

    def workflow_setup(self, name: str, outerdirectory: str = "", **kwargs):
        setupargs, genargs, wargs = filterargs(self, **kwargs)
        self.pre_setup_jobs(outerdirectory=outerdirectory, **setupargs)
        self.setup_jobs(name=name, **genargs, **wargs)

    def workflow_initialize(self, autofetch: bool = True):
        self._auto_fetch(value=autofetch)
        self._initialize_jobs()

    def workflow_hpc_setup(self,
                           server: str,
                           remove: bool = True,
                           backup: bool = False,
                           save_db: bool = True,
                           save_calc: bool = False, **kwargs):

        self.hpcjob_setup(server=server,
                          remove=remove,
                          backup=backup,
                          save_db=save_db,
                          save_calc=save_calc,
                          **kwargs)
        self.submit()

    def hpcjob_setup(self,
                     server: str,
                     remove: bool = True,
                     backup: bool = False,
                     save_db: bool = True,
                     save_calc: bool = False,
                     **kwargs):

        for j in self.jobs:
            hpc = HPCjob(server=server, jobname=j.name, working_directory=j.working_directory)
            j.hpc = hpc
            submitargs, prepargs = filterargs(hpc.server, **kwargs)

            j.submit_hpc(remove=remove,
                         backup=backup,
                         save_db=save_db,
                         dry_run=True,
                         save_calc=save_calc, **submitargs)

            hpc.prepare(**prepargs)

    def submit(self):
        for j in self.jobs:
            j.hpc.submit()

    @staticmethod
    def read_Jobs(wrkpaths:list):
        jobs = []
        for p in wrkpaths:
            j = Gpaw.read_job(working_directory=p)
            jobs.append(j)

        return jobs

    def get_back_jobs(self):
        return self.read_Jobs(wrkpaths=self.wrkpaths)


def filterargs_old(wf: Workflow, **kwargs):

    gen = wf.pre_setup_jobs.__code__
    gen = list(gen.co_varnames[:gen.co_argcount])

    setup = wf.setup_jobs.__code__
    setup = list(setup.co_varnames[:setup.co_argcount])

    setupargs = {}
    genargs = {}
    kwarcopy = {}

    for k, v in kwargs.items():

        if k in gen:
            genargs[k] = v
        elif k in setup:
            setupargs[k] = v
        else:
            kwarcopy[k] = v

    return setupargs, genargs, kwarcopy
