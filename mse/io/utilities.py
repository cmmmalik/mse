
def read_atoms(file:str,
                verbosity:int=1,
                ):
    from ase.io import read

    atoms = read(file, index=":")
    if verbosity >= 1:
        print(atoms[0])
        print("Total Number of atoms read:{}".format(len(atoms)))

    return atoms


def write_to_db(atoms, outdb, **extra_keys):
    from ase.db import connect
    outdb = connect(outdb)
    with outdb:
        for at in atoms:
            outdb.write(at, **extra_keys)