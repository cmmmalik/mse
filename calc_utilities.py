from ase.cell import Cell as asecell
import numpy as np


def get_kpoints_from_density(density: int or float,
                             cell,
                             even: bool = False):

    if not isinstance(cell, asecell):
        cell = asecell(cell)

    r_cell = cell.reciprocal()
    r_cell = np.linalg.norm(r_cell, axis=-1)
    r_cell = 2*np.pi*r_cell
    n_size = r_cell*density
    if even:
        n_size = np.ceil(n_size/2.0)*2.0
    else:
        n_size = np.ceil(n_size)

    return n_size

