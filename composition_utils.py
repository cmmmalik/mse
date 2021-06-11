from preprocessing.atoms import Composition


class MAXcomp:

    def __init__(self, formula:str):

        self._comp = None
        self.formula = formula
        if not self.check_max:
            raise ValueError("Expected MAX phase composition, but got non-MAX composition".format(self.comp.iupac_formula))

    @property
    def formula(self):
        return self._formula

    @formula.setter
    def formula(self, formula):
        self._formula = formula
        self._comp = Composition(formula)

    @property
    def maxelements(self):
        return self.comp.get_mappings()

    @property
    def reduced_comp(self):
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
