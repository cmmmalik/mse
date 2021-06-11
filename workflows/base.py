from ase.atoms import Atoms
from ase.spacegroup.symmetrize import check_symmetry, refine_symmetry
import os

from Simtools.calculator.gpaw_calc import Gpaw
from Simtools.calculator.vasp_calc import VASP
from Simtools.Jobs.job import HPCjob


class Baseworkflow:

    def __init__(self,
                 atoms: Atoms,
                 working_directory: str = None,
                 calculator_type: "gpaw" or "vasp" = "gpaw",
                 dircheck: bool = True,
                 verbosity: int = 1):

        if dircheck:
            if os.path.exists(working_directory):
                raise ValueError("Working directory {} already exists"
                                 "Set 'check = False' to ignore this".format(working_directory))

        self.atoms = atoms
        self.verbosity = verbosity
        if calculator_type not in ["gpaw", "vasp"]:
            raise ValueError("'calculator_type' must set to 'gpaw' or 'vasp', instead of {}".format(calculator_type))
        self._calculator_type = calculator_type

        self._job = None
        self._working_directory = working_directory

    @property
    def atoms(self):
        return self._atoms

    @property
    def working_directory(self):
        return self._working_directory

    @property
    def job(self):
        return self._job

    @atoms.setter
    def atoms(self, at: Atoms):
        if not isinstance(at, Atoms):
            raise ValueError("Given at is an instance of '{}', instead of {}".format(type(at), Atoms))
        self._atoms = at

    @property
    def calculator_type(self):
        return self._calculator_type

    def check_symmetry(self, symprec: float or int = 1e-06):
        return check_symmetry(self.atoms, symprec=symprec)

    def refine_symmetry(self, symprec: float or int = 0.01):
        return refine_symmetry(self.atoms, symprec=symprec)

    def initialize_job(self,
                       name: str,
                       run_type: str = "static",
                       modeinps: dict = None,
                       calcinps: dict = None,
                       **kwargs):

        if not self.working_directory:
            self._working_directory = kwargs.pop("working_directory")

        if self.calculator_type == "gpaw":
            job = Gpaw(name=name,
                       atoms=self.atoms,
                       working_directory=self.working_directory,
                       run_type=run_type,
                       modeinps=modeinps,
                       calcinps=calcinps,
                       **kwargs)

        elif self.calculator_type == "vasp":
            job = VASP(name=name,
                       atoms=self.atoms,
                       working_directory=self.working_directory,
                       run_type=run_type,
                       calcinps=calcinps,
                       )

        self._job = job
        return job

    @job.setter
    def job(self, value):
        if not isinstance(value, Gpaw):
            raise TypeError("Expected an instance of {}, but received {}".format(Gpaw, type(value)))
        self._job = value

    def _auto_fetch(self, value: bool = True):
        self.job.autofetch = value

    def make_ready(self, autofetch: bool = True):
        self._auto_fetch(autofetch)
        self.job.initialize()

    def hpc_setup(self,
                  server: str,
                  remove: bool = True,
                  backup: bool = False,
                  save_db: bool = True,
                  save_calc: bool = False,
                  **kwargs):

        hpc = HPCjob(server=server,
                     jobname=self.job.name,
                     working_directory=self.job.working_directory,
                     )

        self.job.hpc = hpc
        random_folder = kwargs.get("random_folder", True)
        folder_name = kwargs.get("folder_name", False)
        self.job.submit_hpc(remove=remove,
                            backup=backup,
                            save_db=save_db,
                            dry_run=True,
                            save_calc=save_calc,
                            **kwargs)

        hpc.prepare(random_folder=random_folder, folder_name=folder_name)

    def submit(self):
        self.job.submit()

    def hpc_ready_submit(self,
                         server: str,
                         remove: bool = True,
                         backup: bool = False,
                         save_db: bool = True,
                         save_calc: bool = False,
                         **kwargs
                         ):

        self.hpc_setup(server=server,
                       remove=remove,
                       backup=backup,
                       save_db=save_db,
                       save_calc=save_calc,
                       **kwargs)
        self.submit()

    def _show_attr(self, attrs):
        st = "{}".format(getattr(self.job, attrs))
        st = "{}:\n{}".format(attrs, st)
        return st

    def show_inputs(self):
        st = self._show_attr("inputs")
        print(st)
