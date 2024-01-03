from preprocessing.atoms import Aatoms


def rows_strcs(Rows: dict, attach_calculator=False):
    all_strcs = {}
    for k, rows in Rows.items():
        if isinstance(rows, (list or tuple)):
            strcs = [Aatoms(r.toatoms(attach_calculator)) for r in rows]
        else:
            strcs = Aatoms(rows.toatoms(attach_calculator))
        all_strcs[k] = strcs

    return all_strcs


def get_calc_parameters(row):
    calculator_parameters = row.calculator_parameters.copy() #new parameters
    if not calculator_parameters:
        calculator_parameters = row.data.copy() # old parameters--can cause bug
        print("From data of row\n{}".format(calculator_parameters))

    return calculator_parameters


def get_kpts(calculator_parameters, verbosity: int = 0):
    try:
        kpts = calculator_parameters["kpden"]
        kden = calculator_parameters["kden"]
        kpts = {k: kden if v == "kden" else True if v in ['True', True] else False if v in ['False', False] else v for k, v in kpts.items()}

    except KeyError:
        kpts = calculator_parameters["kpts"]

    if verbosity >= 1:
        print("Kpts\n{}".format(kpts))

    return kpts


def get_args(calculator_parameters):
    modeinps = {}
    calcinps = {}
    for k, v in calculator_parameters.items():
        if k == "xc":
            calcinps["xc"] = v
        elif k == "encut":
            modeinps["encut"] = v
    calcinps["kpts"] = get_kpts(calculator_parameters)
    return calcinps, modeinps


def iter_kdens(dbrows: dict):
    all_kdens = {}
    for k, rows in dbrows.items():
        kd = [r.get("kden") for r in rows]
        all_kdens[k] = kd
    return all_kdens