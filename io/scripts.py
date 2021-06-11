import re

#outputs .py files, which can guide calculations
#depracated


class InterfaceGPAW:

    @property
    def calctype(self):

        return self._calctype

    @calctype.setter
    def calctype(self, string):

        if not string in self.avail_types:
            raise ValueError("{} type calculation is not currently supported".format(string))

        self._calctype = string

    @staticmethod
    def _e_volume(posfilename):

        script = """from process_tools import do_e_volume
from analysis_tools import e_volume, get_targetvolume_factor
from ase.parallel import world
from ase.io import read,,
do_e_volume(filename=\"%s\")
,,
if world.rank == 0:,,
    orgatoms = read(\"%s\", format=\"vasp\")
    fname = orgatoms.get_chemical_formula(empirical=True)
    vmin, emin, _ = e_volume(trajfile=\"%%s.traj\" %%fname, outfile=\"%%s\" %%fname)
    get_targetvolume_factor(orgatoms, vmin, outfile=\"%%s.dat\" %%fname)
,,
else:
    pass,,
    ,,""" % (posfilename, posfilename)

        script = script.split(",,")
        script = "\n".join(script)

        return script

    @staticmethod
    def write(outfile, script):

        with open(outfile, "w") as f:
            f.writelines(script)

        print(f"I m finished generating the file {outfile}")

    def e_volume(self, outfile="perform_e_volume.py", script=None):

        if not script:
            script = self._e_volume(self.posfilename)

        self.write(outfile, script)

    @staticmethod
    def _get_convergence_script(filename="BPOSCAR", volumechange="[1.01, 1.02]"):

        script = """from ase.io import read, write, Trajectory
from gpaw import GPAW, PW


#initial parameters
decimtol = 8
verbose = True


#optimization
kden =  3
kpden = {\"density\": kden}
encut = 500
eigensolv = \'rmm-diis\'
filename = \"BPOSCAR\"

atoms = read(filename, format=\"vasp\")

traj = Trajectory(\"volume-change.traj\", \"w\")
calc = GPAW(mode=PW(encut), kpts=kpden, eigensolver=eigensolv, txt=\"static.txt\" )

orgcell = atoms.get_cell(complete=True)

for v in %s:

        newatoms = atoms.copy()
        newatoms.set_calculator(calc)
        newatoms.set_cell(orgcell*v, scale_atoms=True)
        newatoms.get_potential_energy()
        traj.write(newatoms)
    """ % volumechange

        script = script.split(",,")
        script = "\n".join(script)
        return script

    @classmethod
    def static_script_bulk(cls,encut=500, kden=4, outfile="input.py", **inputs):

        params = {"poscarname":"POSCAR",
                  "encut": encut,
                  "kden": kden,
                  "gamma": True,
                  "eigensolver": "rmm-diis",
                  "kpden": {"density": "kden"},
                  "swidth": 0.04,
                  "occuname":"fermi-dirac",
                  "xc":"PBE",
                  "outputtxt":"static.txt"}

        for k,v in inputs.items():
            if k in params:
                params[k] = v
            else:
                raise ValueError(f"Unexpected Argument {k}, arguments can be {params.keys()}")

        if params["gamma"]:
            params["kpden"]["gamma"] = True

        # TODO: This is a workaround, but needs proper solution, for backward compatibility
        try:

            kpdens = kpden_parse_back(params["kpden"]) # doesnot work
            kpdens = "kpden = "+ kpdens
        except KeyError:
            kpdens = "" #setting to empty string

        for k,v in params.items():
            if isinstance(v, str):
                params[k] = "\""+v+"\""

        script = """#created by {_classname}
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

filename = {poscarname}
encut = {encut}
kden = {kden}\n""" + kpdens+  "\n" + """eigensolver = {eigensolver}
xc={xc}
occupations = {{"name":{occuname}, "width": {swidth}}}

atoms = read(filename)

calc = create_calc(txt={outputtxt})
atoms.calc = calc
atoms.get_potential_energy()
"""
        script = script.format(_classname=cls.__class__.__name__, **params)
        cls.write(outfile, script)

    @staticmethod
    def _relaxation_script_slab(structurefilename="POSCAR", **kwargs):

        defaults = {"filename": '"%s"' % structurefilename,
                    "encut": 500,
                    "kpts": {"size": (8, 8, 1), "gamma": True},
                    "eigensolver": "'rmm-diis'",
                    "fmax": 0.001,
                    "relal": "'BFGS'",
                    "parallel": {"sl_auto": "False"},
                    "poissonsolver": {"dipolelyer": "xy"},
                    "swidth": 0.03,
                    "occupations": "{'name': 'fermi-dirac', 'width':swidth}",
                    "mask": [1, 1, 0, 0, 0, 1]}

        for k, v in kwargs.items():
            defaults[k] = v
        scalcparameters = []
        for k, v in defaults.items():
            scalcparameters.append("%s = %s" % (k, v))

        scalcparameters = "\n".join(scalcparameters)
        print("Calculation parameters")
        print(scalcparameters)
        script = """from ase.io import read, write
import gpaw_functions_dev as gf
from gpaw import  GPAW, PW, Mixer, Davidson, FermiDirac 
import numpy as np
from ase.parallel import parprint
from ase.constraints import FixAtoms
from ase.build import surface
from ase.parallel import world
,,
def create_calc(txt="-", dedecut=None):

    calc = GPAW(mode=PW(ecut=encut,
                dedecut=dedecut),
        kpts=kpts,
        xc=xc,
        eigensolver=eigensolver,
        occupations=occupations,
        poissonsolver=poissonsolver,
        parallel=parallel,
        txt=txt)
,,
    return calc
,,
#gpaw parameters\n""" + scalcparameters + """
,,
#reading
if "POSCAR" in filename or ".vasp" in filename:
    atoms = read(filename, format="vasp")

else:

    atoms = read(filename)

atoms.pbc = [True, True, False]
,,
#Ionic relaxation
calc = create_calc(txt="rel-ions.txt")
atoms.set_calculator(calc)
gf.optimize(atoms, reltype="ions", 
        relaxalgorithm=relal,
        fmax=fmax)
,,
#dedecut correction
dedecut = gf.get_dedecut(atoms.copy(),
            encut=encut,
            calc=create_calc(txt="-")
# cell-relaxation
calc = create_calc(txt="rel-cell.txt", dedecut=dedecut)
mxene.set_calculator(calc)
gf.optimize(atoms, reltype="cell",
        relaxalgorithm=relal,
        fmax=fmax,
        mask=mask)
,,
# full relaxation
dedecut = gf.get_dedecut(atoms.copy(),
            encut=encut,
            calc=create_calc(txt="-"))
,,
calc = create_calc(txt="rel-all.txt", dedecut=dedecut)
mxene.set_calculator(calc)
gf.optimize(atoms, reltype="full",
        relaxalgorithm=relal,
        fmax=fmax,
        mask=mask)
,,
e0 = atoms.get_potential_energy()
counter = 0
,,
while True:
,,	
    parprint("Iteration No %s" % counter, flush=True)
    counter += 1
,,
    dedecut = gf.get_dedecut(atoms.copy(),
                encut=encut,
                calc=create_calc(txt="-"))
    calc=create_calc(txt="rel-all-%s.txt" %counter, dedecut=dedecut)
    gf.optimize(atoms, reltype="ful", 
            relaxalgorithm=relal,
            fmax=fmax,
            mask=mask)
,,
    if abs(atoms.get_potential_energy() - e0) < 1e-05:
        parprint("Reached iteration criteria", flush=True)
        parprint("Stopping", flush=True)
        break
    else:
        e0 = atoms.get_potential_energy()
,,
parprint("Writing final relaxed(relaxed_atoms.vasp) atoms", flush=True)
write("relaxed_atoms.vasp", atoms, format="vasp", vasp5=True )"""
        script = script.split(",,")
        script = "\n".join(script)
        return script

    def relaxation_slab(self, outfile="optimize_all.py"):

        script = self._relaxation_script_slab(self.posfilename)
        self.write(outfile, script)

    def convergence_script(self, outfile="conv.py", volumechange="[1.01, 1.02]"):

        script = self._get_convergence_script(self.posfilename, volumechange)
        self.write(outfile, script)

class Pygpaw:
    #Under_developement
    def __init__(self, parameters=None, calcparams=None):

        defaultscalc = {"mode": "PW",
                        "kpts": "kpden",
                        "eigensolver": "eigensolver",
                        "occupations": "occupations",
                        "txt": "txt"}

        inputs_params = {"filename": "poscar"}

        defaultscalc.update(calcparams)
        self.calcparams = defaultscalc

        self.inputs = parameters
        self._firstline = "#Created by {_classname}"
        self._indent = " "*4
        self._nline = "\n"
        self._body1 = self._indent
        self._body2 = self._indent*2
        self._body3 = self._indent*3

    def _imports(self, *fromimportargs):
        defaults = ["from ase.io import read",
                    "from gpaw import GPAW, {mode}".format(mode=self.calcparams["mode"])]
        defaults += fromimportargs
        return defaults

    def def_calc(self, **gpawcalc):


        defs = ["def create_calc(txt=\"-\"):"]
        body_1 = ["calc = GPAW"]



def kpden_parse_back(a:dict): # It is for a work around

    a = str(a)
    a = re.sub(r"\'(kden)\'", r'\1', a)
    return "{"+a+"}"