import re
import warnings

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

from mse.utilities import replace_1_string
import re


class EnhancedComposition(Composition):
    
   

    def __mul__(self, other):
        obj = super(EnhancedComposition, self).__mul__(other)
        return EnhancedComposition(obj)

    @property
    def iupac_formula(self):
        formula = super().iupac_formula
        formula = formula.replace(" ", "")
        return formula

    @property
    def refined_iupac_formula(self):

        formula = self.iupac_formula
        formula = replace_1_string(formula)
        return formula

    def get_withoutel_iupac_formula(self, el):

        formula = self.iupac_formula
        formula = re.sub(el + r"\d*\.?\d*", "", formula)
        return formula

    def get_Ael(self):

        """Get the A atom"""

        try:

            return self._Ael

        except AttributeError:

            groupA = [Element.from_row_and_group(i, 13).symbol for i in range(2, 7)]

            els = [i.symbol for i in self.keys() if i.symbol in groupA]

            if len(els) > 1:

                warnings.warn("More than one A element was found")
            else:
                els = els[0]

            self.Ael = els

            return els

    @property
    def Ael(self):
        return self._Ael

    @Ael.setter
    def Ael(self, value):
        self._Ael = value

    def select_elements_type(self, selection):

        if selection == "transition_metal" or selection == "M":

            return [i.symbol for i in self.keys() if i.is_transition_metal or i.is_lanthanoid or i.is_actinoid]

        elif selection == "GroupA" or selection == "A":

            return [i.symbol for i in self.keys() if i.is_post_transition_metal or i.group in range(13, 17) and i.row > 2]

        elif selection == "X":

            return [i.symbol for i in self.keys() if i.symbol in ["C", "N"] ]#or i.row == 2 and i.group > 12]

        elif selection == "T":

            return [i.symbol for i in self.keys() if i.symbol in ["F", "Cl", "H", "OH", "O"]]

        else:

            raise ValueError("Unknown selection keyword")

    def get_mappings(self, ttype=("M", "A", "X")):
        """ returns dictionary containing mappings of M,A,X as keys to individual elements in the composition """
        return {i: ' '.join(self.select_elements_type(selection=i)) for i in ttype}

    @staticmethod
    def remove_brackets(st):
        from collections import Counter
        dict = Counter(st)
        if dict["("] == dict[")"]:
            pass
        else:
            raise ValueError("{} contains unequal number of '(' and ')'".format(st))
        if not dict["("] == 1:
            raise ValueError(" cannot handle nested brackets")

        b = st.find("(")
        e = st.find(")")
        return st[b + 1, e]


class MAXcomp:

    def __init__(self, formula:str):

        self._comp = None
        self.formula = formula
        if not self.check_max():
            raise ValueError("Expected MAX phase composition, but got non-MAX composition".format(self.comp.iupac_formula))

    @property
    def formula(self):
        return self._formula

    @formula.setter
    def formula(self, formula):
        self._formula = formula
        self._comp = EnhancedComposition(formula)

    @property
    def maxelements(self):
        return self.comp.get_mappings()

    @property
    def reduced_comp(self):
        return self.comp.reduced_composition

    @property
    def reduced_composition(self):
        return self.comp.reduced_composition

    @property
    def comp(self):
        return self._comp

    def check_max(self):
        mapp = self.maxelements
        for k in mapp.keys():
            if not mapp[k]:
                return False
        else:
            try:
                self.get_n()
            except AssertionError:
                return False

        return True

    def get_n(self):
        prim_comp = self.reduced_comp
        mapp = self.maxelements
        quan = prim_comp.get_el_amt_dict()
        n = quan[mapp["X"]]
        assert n + 1 == quan[mapp["M"]]

        return quan[mapp["X"]]


class MXene(MAXcomp):

    def __init__(self, formula):
        try:
            super(MXene, self).__init__(formula=formula)
        except ValueError:
            pass

        if not self.check_mxene():
            raise ValueError("Expected MXene composition: got a non-MXene composition: {}".format(self.comp.iupac_formula.replace(" ", "")))

    def check_mxene(self):
        mapp = self.comp.get_mappings(ttype=("M", "A", "X", "T"))
        if mapp["A"]:
            return False

        for k in mapp.keys():
            if k == "A":
                continue
            if not mapp[k] and k != "T":
                return False
        try:
            self.get_n()

        except AssertionError:
            return False

        return True

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.comp)

    def __getattr__(self, item):
        if hasattr(self.comp, item):
            return getattr(self.comp, item)

        else:
            raise AttributeError(f"{self.__class__.__name__} object as no attribute {item}")


class UnconvMAX(MAXcomp):

    def __init__(self, formula:str, Aelement:str=None): # need to enter the A element
        if not Aelement:
            els = re.findall("[A-Z][a-z]?", formula)
            assert len(els) == 3
            Aelement = els[1]
        self._A = Aelement

        try:
            super(UnconvMAX, self).__init__(formula=formula)
        except ValueError:
            pass



    @property
    def maxelements(self):
        maxelements = super(UnconvMAX, self).maxelements
        if  not maxelements["A"]:
            maxelements["A"] = self._A
        maxelements["M"] = " ".join([i for i in maxelements["M"].split() if i != self._A])

        return maxelements


def generate_MAX_formula(M: str = "Ti",
                         A: str = "Al",
                         X: str = "C",
                         n: list = (1, 2, 3)):
    # build up formula
    formulas = []
    for i in n:
        f = "{M}{i_1}{A}{X}{i}".format(M=M, A=A, X=X, i_1=i + 1, i=i)
        formulas.append(f)
    return formulas


def generate_MAX_formulas(M : list or tuple = ["Ti"],
                          A: list or tuple = ["Al"],
                          X: list or tuple = ["C"],
                          n: list or tuple = (1, 2, 3)):
    formulas = {}
    formulas = {"{}-{}-{}".format(m, a, x): generate_MAX_formula(M=m, A=a, X=x, n=n) for m in M for a in A for x in X}
    return formulas


def get_elements_periodic_table(group: list, row: list):
    R = [Element.from_row_and_group(row=r, group=g).symbol for r in row for g in group]
    return R