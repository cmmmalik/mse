import ase.calculators.calculator as asecalculator
from ase.units import Ha
from ase.parallel import parprint
from importlib import import_module
from copy import deepcopy
import warnings


def getnewcalc(calc: asecalculator, txt="-"):
    try:
        parameters = deepcopy(calc.parameters)
    except Exception as e:
        parprint("Run into an Error:{}".format(e), flush=True)
        parprint("Calc parameters:\n{}".format(calc.parameters), flush=True)
     #   parameters = {k: deepcopy(w) if isinstance(w, (list, tuple, dict, set)) else w for k, w in calc.parameters.items()}
        parameters = calc.parameters.copy()
    calcmodule = import_module(calc.__module__)
    Calc = getattr(calcmodule, calc.__class__.__name__)
    parameters.pop("txt", None)
    return Calc(txt=txt, **parameters)


def setmodecalc(calc: asecalculator, dedecut: str or float or int = None, **kwargs):

    mode = calc.todict()["mode"] # could be expensive, TODO: May be an improvement would be to use 'PW' instance
    orgmode = mode.copy()

    if dedecut is not None:
        mode["dedecut"] = dedecut

    mode.update(kwargs)
    if mode == orgmode: # do nothing if mode is not changed
        warnings.warn("Original and updated mode were not different")

    calc.set(mode=mode) # update the mode


def setmodecalc_v1(calc: asecalculator, dedecut: str or float or int = None, **kwargs):

        mode = calc.parameters["mode"]
        org_mode_dict = mode.copy() if isinstance(mode, dict) else mode.todict()

        if dedecut is not None:
            mode.dedecut = dedecut

        for k, value in kwargs.items():
            if k == "ecut":
                value = value / Ha
            setattr(mode, k, value)

        if mode.todict() == org_mode_dict:  # do nothing if mode is not changed
            warnings.warn("Original and updated mode were not different")