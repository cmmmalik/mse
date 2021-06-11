import re


def sort_formula_string(formula: str, verbosity: int = 0):
    if verbosity >= 1:
        print("Received: {}", format(formula))
    formula_list = re.findall("[A-Z][a-z]?[0-9]*\.?[0-9]*", formula)
    formula_list.sort()
    formula_sorted = "".join(formula_list)
    if verbosity >= 1:
        print("Sorted formula: {}".format(formula_sorted))

    return formula_sorted
