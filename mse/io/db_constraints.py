import warnings


def remove_not_todict(atoms):
    """
    Remove constraints that do not have todict method.
    :param atoms: ase.Atoms
    :return: None
    """
    rmindex = []
    for i, const in enumerate(atoms.constraints):
        if not hasattr(const, "todict"):
            rmindex.append(i)

    # at the moment ase does not have allow FIxsymmetry constraint atoms, storage in ase database,
    # we we are removing it.

    if rmindex:
        warnings.warn("Removing the constraint without todict method, before saving into the database", RuntimeWarning)
        [atoms.constraints.pop(i) for i in rmindex]