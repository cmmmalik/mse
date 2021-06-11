#created by type
# static calculation
from ase.io import read
from gpaw import GPAW, PW

def create_calc(txt="-"):
    calc=GPAW(mode=PW(ecut=encut),
            kpts=kpden,
            xc=xc,
            eigensolver=eigensolver,
            occupations=occupations,
            txt=txt,
            )  
    return calc       

filename = "POSCAR"
encut = 500
kden = 4
kpden = {'density': kden, 'gamma': True}
eigensolver = "rmm-diis"
xc="PBE"
occupations = {"name":"fermi-dirac", "width": 0.04}

atoms = read(filename)

calc = create_calc(txt="static.txt")
atoms.calc = calc
atoms.get_potential_energy()
