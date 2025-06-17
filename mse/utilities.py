import re

"This file contains several funtions usable for string processing..."


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

def replace_1_string_re(s):
    import re
    return  re.sub(r'(?<=[A-Za-z])1(?![\d.])', '', s)


def print_nice(dct: dict, separater="\n", keysep=":", format_spaces=[11,7], nestsep="_"):
    """ returns formatted string representation of nested dictionary, tuples, lists, for printing nicely
     in print function.

    :param dct: dict,tuple, list, that needs to be formatted, may contain nested dict, tuples, lists...
    :param separater: separater specifier, default '\n', that separates entries of dct in formatted output string.
    :param keysep: sepecifier, default ':', that separates key,value pairs of dictionary. It is obviously usable
    for dict only
    :param format_spaces: The width of field for key,value respectively. for lists, tuples etc. only first index is used.
    :param nestsep: str, default _, separates nested list items.
    :return: formatted string.
    """
    # have to use a convertion or boolean {True!s}
    sstr = []
    if isinstance(dct, dict):
        sstr += [f"{k:{format_spaces[0]}}{keysep}{v!s:{format_spaces[1]}}" if not isinstance(v, (
        dict, list, tuple)) else "{0:{3}}{2}{1!s:{4}}".format(k, print_nice(v, separater=",", keysep="=",
                                                                          format_spaces=[1, 1]), keysep, *format_spaces)
                 for k, v in dct.items()]
    elif isinstance(dct, (list, tuple)):
        # sstr += f"{separater}".join(dct)
        sstr += ["{0!s:{1}}".format(i, format_spaces[0]) if not isinstance(i, (dict, list, tuple)) else "{0!s:{1}}".format(
            print_nice(i, separater=nestsep, format_spaces=[1, 1]), format_spaces[0]) for i in dct]
    else:  # it will not be build-in iterable
        sstr += [f"{dct!s:{format_spaces[0]}}"]

    sstr = f"{separater}".join(sstr)
    return sstr