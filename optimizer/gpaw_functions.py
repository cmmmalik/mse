from ase.atoms import Atoms as aseatoms
from ase.io import read
from ase.constraints import UnitCellFilter, StrainFilter
from ase.optimize.bfgs import BFGS
from ase.optimize import QuasiNewton
from ase.optimize.bfgslinesearch import BFGSLineSearch as BFGSLS
from ase.optimize.sciopt import SciPyFminCG as CG, SciPyFminBFGS as ScBFGS
from ase.parallel import parprint
from functools import partial
from os import path
from gpaw import GPAW, PW, FermiDirac

from mse.calculator.ase import getnewcalc, setmodecalc

parprint = partial(parprint, flush=True)


def optimize(atoms,
             reltype: str = 'ions',
             fmax=0.01,
             steps=None,
             relaxalgorithm="BFGS",
             mask=None,
             verbose=True,
             **opargs):
    """
        wrapper function for relaxation

    :param atoms: ase atom object
    :param relax: string (cell, full, "") , type of relaxation,
    :param fmax: number, force criteria
    :param relaxalgorithm: relax algorithm
    :param verbose: bool, default True
    :return: atoms object
    """

    optimizer_algorithms = {"QuasiNewton": QuasiNewton, "BFGS": BFGS, "CG": CG, "ScBFGS": ScBFGS, "BFGSLS": BFGSLS}
    # relaxation algorithms

    if relaxalgorithm in optimizer_algorithms:
        pass
    else:
        raise KeyError("The {0} is invalid or  not found."
                       "The available algorithms are: {1}".format(relaxalgorithm, optimizer_algorithms.values()))

    parprint("Relaxation type: {}\nMask : {}\nOther relaxation parameters {}".format(reltype, mask, opargs))
    if reltype == 'full':

        uf = UnitCellFilter(atoms, mask=mask)
        relax = optimizer_algorithms[relaxalgorithm](uf, logfile="rel-all.log", **opargs)
        if verbose:
            parprint("Full relaxation")

    elif reltype == 'cell':
        sf = StrainFilter(atoms, mask=mask)
        relax = optimizer_algorithms[relaxalgorithm](sf, logfile="rel-cell.log", **opargs)
        if verbose:
            parprint("Cell relaxation only")

    elif reltype == 'ions':  # ionic_relaxation
        relax = optimizer_algorithms[relaxalgorithm](atoms, logfile="rel-ionic.log", **opargs)
        if verbose:
            parprint("Ions relaxation only")

    else:
        raise ValueError("'{}' type relaxation is invalid".format(reltype))

    converged = relax.run(fmax=fmax, steps=steps)
    return converged


def get_modify_calculator(calc = None,
                          encut = None,
                          dedecut = None,
                          kp = None,
                          resetpw = False,
                          txt = "-",
                          **kwargs):

    """

    :param encut: energy cutoff
    :param kp: Dictionary i.e. {"density": No., "gamma" : True, "even": True }
    :param gammacentered: bool,
    :param eigensolver: string, type of sigensolver for Hamiltonian
    :param swidth: int, smearing width
    :param psp: Dictionary, to select a different psedupotential, i.e. {"Mg" : "2"}
    :param spin: True, False, spin-polarized/spin-paried calculations
    :param U: Dictionary: +U approach
    :param dedecut:  derivative with of energy respect to energy cut-off  for  pulay correction, number(fload, number)
    :param poissonsolver: Dict, includes - dipole correction for unsymmetrical slab, {"dipolelayer", "xy"}
    :param mixer: tuple, for density mixing (type of mixer, (values))
    :param calc: calculator object
    :param output: output file name (".txt")
    :return: calculator object
    """

    par_args = ( "xc",
                   "swidth",
                   "mixer",
                   "eigensolver",
                   "psp",
                   "spin",
                   "U",
                   "poissonsolver" ,
                   "parallel") # not being used anymore

    pwargs = ("cell",
              "gammacentered",
              "pulay_stress",
              "force_complex_dtype" )


    parameters = {}
    pwparameters = {}

    for key, value in kwargs.items():

        if key in pwargs:
            pwparameters[key] = value
            initializepw = True

        else:
            parameters[key] = value

    if kp is not None:
        parameters["kpts"] = kp

    parameters["txt"] = txt

    if encut is not None or dedecut is not None:
        resetpw = True

    if dedecut:
        assert encut

    if resetpw:
        modePW = PW(encut, **pwparameters)

    if calc is None:
        calc = GPAW(mode=modePW, **parameters)
        return calc

    elif resetpw == True:   # Used Calculator
        parprint("Setting calc mode", flush=True)
        calc.set(mode=modePW)

    # if kp is not None:
    #     calc.set(kpts=kp)
    #     parprint("setting kpts to %s" % kp, flush=True)


    parprint("Setting Parameters:\n{}".format(parameters))
    calc.set(**parameters)



    # if parameters["U"] is not None:
    #     calc.set(setups=parameters["U"])
    #
    # if parameters["psp"] is not None:
    #     calc.set(setups=parameters["psp"])
    #     parprint("setting pseudopotential(psp) %s" %psp, flush=True)
    #
    # if parameters["swidth"] is not None:
    #     calc.set(occupations=FermiDirac(parameters["swidth"]))
    #     parprint("setting smearing-width: %s" % parameters["swidth"])
    #
    # if parameters["spin"] is not None:
    #     calc.set(spinpol=True)
    #     parprint("setting spin to %s" % parameters["spin"], flush=True)
    #
    # if parameters["poissonsolver"] is not None:
    #     calc.set(poissonsolver=parameters["poissonsolver"])
    #
    # if parameters["mixer"] is not None:
    #     attach_mixer(calc, parameters["mixer"])
    #     parprint("calling the mixer function:", flush=True)


def attach_mixer(calc, mixer):

    """To change default values of mixer parameters of calculator
    :parameter:

            calc: Calculator Object
            mixer: nested tuple, ("string", (1,2,3))

    """
    from gpaw import Mixer, MixerDif,  MixerSum, MixerSum2

    mixtypedict = {"Mixer": Mixer,
                   "MixerDif": MixerDif,
                   "MixerSum": MixerSum,
                   "MixerSum2": MixerSum2}

    if mixer[0] in mixtypedict:
        calc.set(mixer=mixtypedict[mixer[0]](**mixer[1]))

    else:

        raise Exception("%s is invalid/incorrect" % mixer[0])


def get_dedecut(atoms,
                encut,
                calc,
                **kwargs):


    """
    :param atoms: ase atom object, required
    :param encut: energy cut-off, required
    :param kp: Dictionary i.e. {"density": No., "gamma" : True, "even": True }
    :param calc: calculato object, required
    :param decimtol: decimal tolerance for rounding the pulay stress , default 8,
    :param verbose: silent or visible ,True, False
    :return: dedecut value for pulay correction
    """
    defaultparameters = { "verbose": True,
                          "suffixf": "",
                          "step": 10}

    for key, value in kwargs.items():
        if key in defaultparameters:
            defaultparameters[key] = value # updating the dictionary
        else:
            raise KeyError("Unkownd/invalid keyword argument %s" % key)

    e = []
    step = defaultparameters["step"]
    verbose = defaultparameters["verbose"]

    for icut in [encut-step, encut+step]:
        ats = atoms.copy()
        parprint("Encut in loop is %s" % icut, flush=True)
        get_modify_calculator(encut=icut,
                              calc=calc,
                              txt="dedecut.txt".format(icut, defaultparameters["suffixf"]))

        ats.set_calculator(calc)
        e.append(ats.get_potential_energy())

    dedecut = (e[1]-e[0])/(2*step)

    if verbose:
        parprint("energy at two different cutoffs: {}".format(e), flush=True)
        parprint("Pulay stress correction is {}".format(dedecut), flush=True)

    return dedecut


def dedecut_outputfile(encut):

    firstname = "deduct"
    counter = 0

    while True:

        fullname = "%s-%s" % (firstname, counter) + "-%s.txt " % encut
        if not path.isfile(fullname):  # checking in the current directory only
            break
        counter += 1
        parprint("dedecut filename is: %s" % fullname, flush=True)

    return fullname


def dedecut_read(fileu, filel, step=10, decimtol=8, verbose=True):

    atomsu = read(fileu)
    atomsl = read(filel)

    dedecut = (atomsu.get_potential_energy() - atomsl.get_potential_energy())/2*step
    dedecut = round(dedecut, decimtol)

    if verbose:
        parprint("The pulay correction using dedecut is %s " % dedecut, flush=True)

    return dedecut


def get_dedecut_v1(atoms: aseatoms ,
                encut: int or float,
                calc = None,
                step_encut: int = 10,
                decimtol: int = 8,
                **kwargs):


    """
    :param atoms: ase atom object, required
    :param encut: energy cut-off, required
    :param kp: Dictionary i.e. {"density": No., "gamma" : True, "even": True }
    :param calc: calculato object, required
    :param decimtol: decimal tolerance for rounding the pulay stress , default 8,
    :param verbose: silent or visible ,True, False
    :return: dedecut value for pulay correction
    """
    defaultparameters = { "verbose": True,
                          "suffixf": "",
                          "prefixf":"dedecut"}

    assert atoms.calc or calc

    for key, value in kwargs.items():
        if key not in defaultparameters:
            raise KeyError("Unkown/invalid keyword argument %s" % key)

    defaultparameters.update(kwargs)
    verbose = defaultparameters["verbose"]

    if not calc:
        calc = getnewcalc(atoms.calc) # do not change the original calculator, creat new calculator

    parprint("Starting dedecut-estimation calculations")
    energy = []

    tailname = "{}".format("-"+defaultparameters["suffixf"])

    atoms.calc = calc

    for icut in [encut-step_encut, encut+step_encut]:
        parprint("Calculation with encut:{}".format(icut))
      #  calc = get_modify_calculator(encut=icut,
      #                               calc=calc,
      #                               txt="dedecut-{}-{}.txt".format(icut, defaultparameters["suffixf"]))

        filename = "dedecut-{}{}.txt".format(icut, tailname)
        setmodecalc(calc=calc, dedecut=None, ecut=icut)
        calc.set(txt=filename)
        calc.reset()
        energy.append(atoms.get_potential_energy())

    dedecut = (energy[1]-energy[0])/(2*step_encut)
    dedecut = round(dedecut, decimtol)

    if verbose:
        parprint("energy at two different cutoffs", energy, flush=True)
        parprint("pulaty correction is %s" % dedecut, flush=True)

    return dedecut
