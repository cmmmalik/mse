from ase.cell import Cell as asecell
import numpy as np

import warnings


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
        n_size = (np.ceil(n_size/2.0)*2.0).astype(int)
    else:
        n_size = np.ceil(n_size).astpye(int)

    return n_size


def get_density_from_kpoints(kpts:list or tuple,
                             cell,
                             even:bool=False):
    if not isinstance(cell, asecell):
        cell = asecell(cell)

    if even:
        kpts = np.floor(kpts/2.0)*2.0

    r_cell = 2*np.pi*cell.reciprocal().lengths()
    density = np.asarray(kpts)/r_cell

    density = np.floor(density).astype(int)

    if len(np.unique(density)) != 1:
        warnings.warn("density values are different....")
        print(density)
        return np.min(density)

    return density[0]


def get_size_ak(n_ak:int or float,
                cell,
                axis:int=0,
                even:bool=None):
    """
    Takes grid size(number of kpoints) along principle axis (given by 'axis', e.g. 0 for a lattice vector)
    and gives appropriate gride size of kpoints by calculating size along other two principle direction.
     .
    :param n_ak: int
    :param cell: atoms cell.
    :return:
    """
    lengths = cell.lengths()
    na = n_ak*lengths[axis]
    size = np.round(na/lengths).astype(int)
    if even is not None:
        if even:
            size = (np.ceil(size/2.0)*2.0).astype(int)
        else:
            size = np.ceil(size).astype(int)
    # size = np.ceil([n_ak, n_bk, n_zk])
    return size