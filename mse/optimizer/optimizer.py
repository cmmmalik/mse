import ase.atoms
from itertools import cycle as itcycle
from typing import Iterable
import warnings

from mse.calculator.ase import getnewcalc, setmodecalc, resetmodecalc
from mse.optimizer.gpaw_functions import get_dedecut as gget_dedecut, optimize as goptimize

import json
from monty.json import MontyDecoder, MontyEncoder

# TODO: test this.
#  expecting changes in the atoms here, should also reflect in the upper levels.
# TODO: Merge gpaw_functions into optimize class or in calculators, remove the dependency from them.
# TODO: Create Documentation for showing the usage and pointing out the problems and
#  additional possible required improvements


class Optimize:
    """
     A front-end class interface of function gpaw_functions.optimize
    """

    _allowedtypes = ["ions", "cell", "full"]
    default_parameters = (("fmax", 0.01,),
                          ("algorithm", "BFGS",),
                          ("trajectory", None,),
                          ("mask", None,),
                          ("db", None))

    def __init__(self, atoms: ase.atoms,
                 reltype: str or None = "pos",
                 working_directory:str=".",
                 extra_algorithmkwargs=None,
                 **kwargs):

        if extra_algorithmkwargs is None:
            extra_algorithmkwargs = {}
        self._atoms = atoms
        self.reltype = reltype
        self.parameters = dict(self.default_parameters)
        self.working_directory=working_directory

        # aargs = self.allowed_args()

        # for k in kwargs:
        #    if k not in aargs:
        #        display = ",".join(aargs)
        #        raise ValueError(f"Invalid keyword argument {k}, allowed arguments are {display}")

        self.steps = kwargs.get("steps", None)
        upde = None
        if self.steps:
            upde = True

        self._update_dedecut = kwargs.get("update_dedecut", upde)
        self.cycle = kwargs.get("cycle", False)

        for k, v in self.parameters.items():
            self.parameters[k] = kwargs.get(k, v)
        if kwargs.get("db", None):
            self.parameters["db"] = "optimization.db"
        self.dedecut_calc = kwargs.get("dedecut_from_calc", False)
        self.dedecut_stepsize = kwargs.get("dedecut_stepsize",  10)

        if not extra_algorithmkwargs:
            extra_algorithmkwargs = dict()
        self._extra_algo_kwargs = extra_algorithmkwargs

    @classmethod
    def allowed_args(cls):
        args, _ = zip(*cls.default_parameters)
        return ["steps", "update_dedecut", "cycle"] + list(args)

    @property
    def atoms(self):
        return self._atoms

    @property
    def cycle(self):
        return self._cycle

    @cycle.setter
    def cycle(self, value):
        if isinstance(value, bool):
            self._cycle = value
        else:
            raise ValueError(f"Expected an instance of {bool}, but got {type(value)}")

    @property
    def steps(self):
        return self._steps

    @steps.setter
    def steps(self, value):
        if isinstance(value, Iterable) and not isinstance(value, str) or value is None:
            self._steps = value
        else:
            raise ValueError(f"Expected an instance of {Iterable}, but got {type(value)}")

    @property
    def reltype(self):
        return self._reltype

    @reltype.setter
    def reltype(self, value):

        if value not in self._allowedtypes:
            raise ValueError(f"Unexpect value of reltype {value}, it  must be one of {self._allowedtypes}")
        self._reltype = value

    def relax(self, **kwargs):

        kwargs.update(self.parameters)
        kwargs.pop("algorithm", None)
        kwargs.update(self._extra_algo_kwargs)
        # because the name is different in the function call. due to backward compatibility reasons
        conv = goptimize(self.atoms,
                         reltype=self.reltype,
                         relaxalgorithm=self.parameters["algorithm"],
                         working_directory=self.working_directory,
                         **kwargs)
        return conv

    def get_dedecut(self):

        ncalc = getnewcalc(self.atoms.calc, txt="dedecut.txt")
        # Only for gpaw this will work
        encut = ncalc.todict()["mode"]["ecut"]
        dedecut = gget_dedecut(self.atoms.copy(), encut=encut, calc=ncalc, step=self.dedecut_stepsize)
        return dedecut

    def optimize(self):

        if self.reltype == "ions": # no need of dedecut calculation in this case

            return self._optimize_pos()

        else:
            return self._optimize()

    def _optimize_pos(self):

        if self.steps:
            generator = self.steps
            if self.cycle:
                generator = itcycle(self.steps)
            return self._optimize_recursive_static(generator=generator)

        return self.relax()

    def _optimize_recursive_static(self, generator):
        for c, st in enumerate(generator):
            conv = self.relax(steps=st)
            if conv:
                return conv

    def _check_apply_dedecut(self):
        calc_dedecut = getattr(self.atoms.calc.parameters["mode"], "dedecut", None)
        if not calc_dedecut:
            dedecut = self.get_dedecut()
            self._updatemodecalc(dedecut=dedecut)

        else:
            if calc_dedecut == "estimate":
                warnings.warn(
                        "The calculator built-in method will be used for getting the dedecut(derivative of energy w.r.t encut)")

            else:
                warnings.warn(f"dedecut was provided as, {calc_dedecut}")

    def _optimize(self):

        if not self.steps: # ignore dedecut_update, do only once before calling the actual optimization function

            self._check_apply_dedecut()
                # dedecut = self.get_dedecut()
                # self._updatemodecalc(dedecut=dedecut)

            conv = self.relax()
            return conv

        if self.cycle: #cycle over the steps over and over again, till convergence is achieved.
            generator = itcycle(self.steps)
        else:
            generator = self.steps
        return self._optimize_recursive(generator=generator)

    def _optimize_recursive(self, generator):

        for c, st in enumerate(generator):

            if self.update_dedecut or c == 0:

                self._check_apply_dedecut()
                # dedecut = self.get_dedecut()
                # self._updatemodecalc(dedecut=dedecut)
            conv = self.relax(steps=st)

            if conv:
                return conv

    def _updatemodecalc(self, dedecut, **kwargs):
        resetmodecalc(self.atoms.calc, dedecut=dedecut, **kwargs)

    @property
    def update_dedecut(self):
        return self._update_dedecut

    @update_dedecut.setter
    def update_dedecut(self, value):
        if isinstance(value, bool):
            self._update_dedecut = value
        else:
            raise ValueError(f"Expected an instance of {bool}, but got {type(value)}")

    def todict(self):
        dct = {}
        for k,v in self.__dict__.items():
            if isinstance(v, ase.atoms.Atoms):
                v=v.todict()
            dct[k] = v

        return dct

    def as_dict(self):
        dct = {"@module": self.__class__.__module__,"@class": self.__class__.__name__}
        dct.update(self.todict())
        return dct

    @classmethod
    def from_dict(cls, dct):

        obj = cls(atoms=ase.atoms.Atoms.fromdict(MontyDecoder().process_decoded(dct["_atoms"])),
                  reltype=dct["_reltype"],
                  working_directory=dct.get("working_directory"),
                  )
        for k,v in dct.items():
            if k in ["_atoms", "_reltype", "working_directory"]:
                continue
            setattr(obj, k, v)
        return obj

    def tojson(self):
        return json.dumps(self, cls=MontyEncoder)

    @classmethod
    def fromjson(cls, jsonstring):
        obj = json.loads(jsonstring, cls=MontyDecoder)
        return obj


class Moptimize:

    def __init__(self, atoms: ase.atoms,
                 working_directory:str=".",
                 reltype: Iterable = ("ions", "cell", "full"), **kwargs):

        self._atoms = atoms
        self.parameters = kwargs
        self.reltype = reltype
        self.working_directory=working_directory

    def optimizers(self):
        for i, r in enumerate(self.reltype):
            kwargs = {}
            for key, value in self.parameters.items():
                if isinstance(value, (list, tuple)):
                    kwargs[key] = value[i]
                else:
                    kwargs[key] = value
            yield Optimize(atoms=self.atoms, reltype=r, working_directory=self.working_directory, **kwargs)

    def optimize(self):
        for op in self.optimizers():
            conv = op.optimize()
            if not conv:
                warnings.warn("The relaxation_type {} did not converge".format(op.reltype), RuntimeWarning)
        return conv

    @property
    def reltype(self):
        return self._reltype

    @reltype.setter
    def reltype(self, value):
        if isinstance(value, str):
            value = [value]
        self._reltype = value

    @property
    def atoms(self):
        return self._atoms

    def to_dict(self):
        dct = {}
        for k,v in self.__dict__.items():
            if isinstance(v, ase.atoms.Atoms):
                v=v.todict()
                if "constraints" in v:
                    const_lst = []
                    for c in v["constraints"]:
                        try:
                            const_lst.append(c.todict())
                        except AttributeError:
                            new_dct = {"name": c.__class__.__name__, # that's whi ase knows which class to import
                                        }
                            new_dct.update(c.__dict__)
                        if new_dct["name"] == "FixSymmetry": # at the moment  atoms are also inserted in the dictionary..
                            new_dct["kwargs"].pop("atoms")
                        const_lst.append(new_dct)

                    v["constraints"] = const_lst
            dct[k] = v

        return dct

    def as_dict(self):
        dct = {"@module": self.__class__.__module__, "@class": self.__class__.__name__}
        dct.update(self.to_dict())
        return dct

    @classmethod
    def from_dict(cls, dct):
        dct = MontyDecoder().process_decoded(dct)
        obj = cls(atoms=ase.atoms.Atoms.fromdict(dct["_atoms"]),
                  working_directory=dct["working_directory"],
                  reltype=dct["_reltype"])

        for k,v in dct.items():
            if k in ["_atoms", "_reltype", "working_directory"]:
                continue

            setattr(obj, k, v)
        return obj

    def tojson(self):
        return json.dumps(self, cls=MontyEncoder)

    @classmethod
    def fromjson(cls, jsonstring):
        obj = json.loads(jsonstring, cls=MontyDecoder)
        return obj




