from ase.atoms import Atoms as aseatoms
from ase.db.row import AtomsRow
from pymatgen.core.composition import Composition as Pymatcomposition
import warnings

from Simtools.database_utilities import Search as dbsearch
from preprocessing.atoms import Aatoms

# Do the analysis(make it ready for the next step)  on or of energies-related quantities


class Analyzer:

    def __init__(self,
                 formula: str,
                 atoms: aseatoms = None,
                 aatoms: Aatoms = None,
                 row: AtomsRow = None,
                 attach_calculator=False):

        assert atoms or aatoms or row

        self._atoms = None
        self._row = None
        self.formula = formula

        if atoms:
            self.atoms = atoms

        if isinstance(row, AtomsRow):
            self._row = row
            if not self.atoms:
                self.atoms = row.toatoms(attach_calculator)

    @property
    def aatoms(self):
        return Aatoms(self._atoms)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, ats: aseatoms):
        if isinstance(ats, aseatoms):
            self._atoms = ats
        else:
            raise ValueError("Expected an instance of {}, but got {}".format(aseatoms, type(aseatoms)))

    @property
    def formula(self):
        return self._formula

    @formula.setter
    def formula(self, value: str):
        if isinstance(value, str):
            self._formula = value
        else:
            raise ValueError("Expected an instance of {}, but got {}".format(str, type(value)))

    @property
    def row(self):
        return self._row

    @row.setter
    def row(self, value: AtomsRow):
        if isinstance(value, AtomsRow):
            self._row = value
            rowatoms = value.toatoms(attach_calculator=False)
            if not self.atoms:
                self.atoms = rowatoms
            elif self.atoms != rowatoms:
                warnings.warn("Atoms object in the attribute row:{} and atoms{} differs".format(rowatoms, self.atoms))
        else:
            raise ValueError("Expected an instance of {}, but got {}".format(AtomsRow, type(value)))

    def get_elements(self):
        els = list(self.aatoms.get_el_amt_dict().keys())
        return els

    def search_in_database(self, db, *args, **kwargs):
        searcher = dbsearch(db=db)
        rows = searcher.get_rows_formula(formula=self.formula, *args, **kwargs)
        return rows

    def search_elements_in_database(self, db, *args, **kwargs):
        searcher = dbsearch(db)
        els = self.get_elements()
        els_rows = {}
        for element in els:
            rows = searcher.get_rows_formula(formula=element, *args, **kwargs)
            if len(rows) >= 1:
                warnings.warn("More than one rows are found for the element {}".format(element))
            els_rows[element] = rows
        return els_rows


def energy_per_atom(row=None,
                       atoms: aseatoms = None):
    assert row or atoms
    if row:
        energy = row.energy
        natoms = row.natoms
    else:
        energy = atoms.get_potential_energy()
        natoms = len(atoms)
    return energy/natoms


def energy_per_formula(row=None,
                       atoms: aseatoms=None):

    en_per_atom = energy_per_atom(row=row, atoms=atoms)
    if row:
        formula = row.formula
    else:
        formula = atoms.get_chemical_formula(empirical=True)
    comp = Pymatcomposition(formula).reduced_composition
    return en_per_atom*comp.num_atoms



