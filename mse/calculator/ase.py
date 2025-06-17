import ase.calculators.calculator as asecalculator
from ase.units import Ha
from gpaw import PW
from ase.parallel import parprint
from importlib import import_module
from copy import deepcopy
import warnings
import sys


def getnewcalc(calc: asecalculator, txt="-"):
    try:
        parameters = deepcopy(calc.parameters)
        parallel = deepcopy(calc.parallel)
    except Exception as e:
        parprint("Run into an Error:{}".format(e), flush=True)
        parprint("Calc parameters:\n{}".format(calc.parameters), flush=True)
     #   parameters = {k: deepcopy(w) if isinstance(w, (list, tuple, dict, set)) else w for k, w in calc.parameters.items()}
        parameters = calc.parameters.copy()
    calcmodule = import_module(calc.__module__)
    Calc = getattr(calcmodule, calc.__class__.__name__)
    parameters.pop("txt", None)
    return Calc(txt=txt, **parameters, parallel=parallel)


def setmodecalc(calc: asecalculator, dedecut: str or float or int = None, **kwargs):

    mode = calc.todict()["mode"] # could be expensive, TODO: May be an improvement would be to use 'PW' instance
    orgmode = mode.copy()

    if dedecut is not None:
        mode["dedecut"] = dedecut

    mode.update(kwargs)
    if mode == orgmode: # do nothing if mode is not changed
        warnings.warn("Original and updated mode were not different")

    calc.set(mode=mode) # update the mode


def resetmodecalc(calc: asecalculator, dedecut: str or float or int = None, **kwargs):

        mode = calc.parameters["mode"]
        org_mode_dict = deepcopy(mode) if isinstance(mode, dict) else mode.todict()


        if "pw" not in org_mode_dict.values():
            warnings.warn("dedecut can not be set for non PW method, I will just update with kwargs")
            org_mode_dict.update(kwargs)
            calc.set(mode=org_mode_dict.copy())
            return

        mode_dct = deepcopy(org_mode_dict)
        mode_dct.update(kwargs)

        mode_dct["dedecut"] = dedecut
        mode_dct.pop("name", None)
        mode = PW(**mode_dct)

        # if dedecut is not None:
        #     mode.dedecut = dedecut
        #
        # for k, value in kwargs.items():
        #     if k == "ecut":
        #         value = value / Ha
        #     setattr(mode, k, value)

        if mode.todict() == org_mode_dict:  # do nothing if mode is not changed
            warnings.warn("Original and updated mode were not different")
        parprint("Resetting the mode to: {}".format(mode.todict()), flush=True)
        calc.set(mode=mode)