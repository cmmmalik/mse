from ase.db import core as dBcore
from ase.atoms import Atoms as aseatoms
from ase.spacegroup.symmetrize import check_symmetry
from collections import Iterable
from pymatgen.io.ase import AseAtomsAdaptor

from database_utilities import get_rows_db
from preprocessing.atoms import Pycomp
from Simtools.analysis.symmetry import get_symmetry_kp_gamma
from Simtools.ext.materials_project import get_entries as matentries
from Simtools.workflows.kpts_new import Workflow as KPworkflow


class Automatized_worker:

    def __init__(self, compound: str, atoms: aseatoms= None, verbosity: int=1):

        self._compound = compound
        self._atoms = atoms
        self.verbosity = verbosity

        self._mat_entry = None
        self._kp_workers = None

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, value: aseatoms):
        if isinstance(value, (aseatoms, None)):
            self._atoms = value

    @property
    def compound(self):
        return self._compound

    @property
    def kp_workers(self):
        return self._kp_workers

    @property
    def mat_entry(self):
        return self._mat_entry

    def in_database(self, db: str or dBcore, *args, **kwargs):
        dbrows = get_rows_db(db, self.compound, *args, **kwargs)
        dbrf = [r for r in dbrows if Pycomp(r.formula).reduced_composition == Pycomp(self.compound).reduced_composition]
        return dbrf

    def check_in_database(self, db: str or dBcore, *args, **kwargs):
        dbrf = self.in_database(db, *args, **kwargs)
        if dbrf:
            return True
        return False

    def materials_project_entry(self, mpkey: str, **kwargs):
        return matentries(formula=self.compound, mpkey=mpkey, **kwargs)

    def set_mat_entry(self, mpkey:str, sort_by_e_above_hull=True, **kwargs):
        self._mat_entry = matentries(formula=self.compound,
                                     mpkey=mpkey,
                                     sort_by_e_above_hull=sort_by_e_above_hull,
                                     **kwargs)

    def get_mat_atoms_entry(self):
        assert self.mat_entry
        atomss = [AseAtomsAdaptor.get_atoms(entry) for entry in self.get_mat_atoms_entry()]
        return atomss

    def kp_convergence(self,
                      encut: int or float,
                      kplimits: list or tuple or Iterable,
                      verbosity: int = 1):

        assert self.atoms
        workers = KPworkflow(encut=encut,
                             atoms=self.atoms,
                             kplimits=kplimits,
                             verbosity=verbosity)
        self._kp_workers = workers
        return workers

    def kp_gamma(self, symprec: float = 1.0e-6, verbose: bool = False):
        symmetry = check_symmetry(atoms=self.atoms, symprec=symprec, verbose=verbose)
        return get_symmetry_kp_gamma(symmetry=symmetry)

    def set_kp(self,
              symmetrize: bool = True,
              autofetch: bool = True,
              **kwargs):

        kp_gamma = self.kp_gamma()
        self.kp_workers.workflow_setup(name=self.compound,
                                       outerdirectory=self.compound,
                                       symmetrize= symmetrize,
                                       kpts=kp_gamma,
                                       **kwargs)
        if self.verbosity >= 1:
            print(self.kp_workers.jobs)

        self.kp_workers.workflow_initialize(autofetch=autofetch)

    def submit_kp(self,
                  server: str,
                  remove: bool = True,
                  backup: bool = False,
                  save_db: bool = True,
                  save_calc: bool = False,
                  **kwargs):

        self.kp_workers.hpcjob_setup(server=server,
                                     remove=remove,
                                     backup=backup,
                                     save_db=save_db,
                                     save_calc=save_calc,
                                     **kwargs)
        self.kp_workers.submit()

    def analyzer(self):
        raise NotImplementedError("Under developement")



