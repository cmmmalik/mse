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


def replace_1_string(s):
    plindex = []

    for i, f, m, l in zip(range(len(s)), s, s[1:], s[2:]):

        if m == "1" and f.isalpha() and l.isalpha():
            plindex.append(i + 1)
    try:
        int(s[-2:])
    except ValueError:  # Means the last integer

        try:
            float(s[-3:])

        except ValueError:

            if s[-1] == "1":
                plindex.append(len(s) - 1)
    strwithoutone = "".join([value for i, value in enumerate(s) if i not in plindex])

    return strwithoutone

