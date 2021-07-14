import os

from ase.db import connect
from ase.db.core import Database as dBcore
from collections import Iterable

from preprocessing.atoms import Composition as Pycomp


def get_rows_db(db: str or dBcore = "output.db", *query, **kwargs):

    def _iter_rows(mydb, *query, **kwargs):
        Rows = []
        for row in mydb.select(*query, **kwargs):
            Rows.append(row)
        return Rows

    if isinstance(db, str):
        with connect(db) as mydb:
            return _iter_rows(mydb, *query, **kwargs)

    return _iter_rows(db, *query, **kwargs)


def formulas(rows):
    return [r.formula for r in rows]


def primitive_formula(rows):
    return [Pycomp(r.formula).reduced_composition for r in rows]


def get_fromdatabase(formula, database):
    if not os.path.exists(database):
        raise FileNotFoundError("The database does not exist")

    with connect(database) as mydb:

        if isinstance(formula, str):
            return mydb.get_atoms(formula=formula)
        elif isinstance(formula, Iterable):
            out = {}
            for f in formula:
                rows = [r for r in mydb.select(f)]
                out[f] = rows

    return out


def iter_getfromdatabase(elements, database):
    out = {}
    for k, v in elements.items():
        rows = get_fromdatabase(formula=v, database=database)
        out[k] = rows

    return out


class Search:

    def __init__(self, db: dBcore or str):
        if not isinstance(db, (dBcore, str)):
            raise ValueError("Expected an instance of {}, but got {}".format(dBcore, type(db)))

        if isinstance(db, str):
            if not os.path.exists(db):
               raise FileNotFoundError("The database does not exist")
        self._db = db

    @property
    def db(self):
        return self._db

    def get_rows(self, *args, **kwargs):
        return get_rows_db(db=self.db, *args, **kwargs)

    def get_rows_formula(self, formula, *args, **kwargs):
        rows = self.in_database(formula, *args, **kwargs)
        #re-check
        rows = [r for r in rows if Pycomp(r.formula).reduced_composition == Pycomp(formula).reduced_composition]
        return rows