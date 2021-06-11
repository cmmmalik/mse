from collections import OrderedDict
import numpy as np
from scipy.linalg import solve, lstsq, LinAlgError, svdvals
from itertools import combinations as itcomb

from preprocessing.atoms import Composition as Pycomp


class LinearlydependentMatrix(Exception):
    pass


def equation_balancer(reactants: list,
                      products: list,
                      verbosity: int = 1,
                      decimal_tolerance: int = 4):
    '''Balances a given chemical equation,
    input parameters:
        reactants, products, verbosity(default=1)'''

    def fill_inmatrix(matrix, compounds: list, elsindex: np.asarray):
        for col, reacs in enumerate(compounds):
            for els in reacs:
                row = (elsindex == els.name).nonzero()
                assert len(row) == 1 and len(row[0]) == 1
                row = row[0][0]
                quan = reacs[els]
                matrix[row][col] = quan
        if verbosity >= 2:
            print("Debug:mode\n{Matrix}")
            print(matrix)
        return matrix

    org_reactants = reactants
    org_products = products

    reactants = [Pycomp(r) if isinstance(r, str) else r for r in reactants]
    products = [Pycomp(p) if isinstance(p, str) else p for p in products]

    lhs_els = np.unique([el.name for r in reactants for el in r])
    lhs_els.sort()
    rhs_els = np.unique([el.name for p in products for el in p])
    rhs_els.sort()

    if (lhs_els != rhs_els).all():
        raise ValueError("Reactants and products do not contain equal type of number of unique  elements")

    n = 1
    a = np.zeros((len(lhs_els), len(reactants)))
    b = np.zeros((len(rhs_els), len(products)))
    # fill_in
    a = fill_inmatrix(matrix=a, compounds=reactants, elsindex=lhs_els)
    b = fill_inmatrix(matrix=b, compounds=products, elsindex=rhs_els)
    orgb = b.copy()
    orga = a.copy()

    ab = -1 * b[:, 1:]

    try:
        a = np.hstack((a, ab))
    except ValueError:
        a = np.hstack((a, ab[np.newaxis].T))

    b = b[:, 0]
    assert b is not None

    if verbosity >= 2:
        print("Debug:\nSolving the MAtrix:\n{}X={}".format(a, b))

    # Now Solve
    try:
        coeffs = solve(a, b)
    except ValueError:
        rows, columns = a.shape
        if rows >= columns:
            diff = columns - rows
            coeffs = solve(a[:diff], b[:diff])
            coeffs2 = solve(a[-diff:], b[-diff:])
            if verbosity >= 2:
                print("Debug:\n{}".format(coeffs))
                print("Coeffs_2:{}".format(coeffs2))
            if not np.allclose(coeffs, coeffs2):
                raise RuntimeError("Can't balance the equation")
    except LinAlgError:  # non-singular Matrix
        coeffs = lstsq(a, b)[0]
        print(coeffs)

    #            try:
    #                coeffs,_,_,_ = lstsq(a,b)
    #                coeffs = np.around(coeffs, 3)
    #           except LinAlgError:
    #               raise LinAlgError("Couldnt balance the equation")
    try:
        factor = np.gcd.reduce(coeffs)
    except Exception:
        factor = 1

    coeffs = coeffs / factor
    n = n / factor
    if verbosity >= 2:
        print("Debug:\n{}".format(np.absolute(coeffs)))

    coeffs = np.around(coeffs, decimal_tolerance)

    if (coeffs == 0).any(-1):
        raise RuntimeError("One of the coefficient is zero: {}".format(coeffs))
    zz = zip(coeffs, reactants + products[1:])
    if verbosity:
        print("Debug:\n{}".format(coeffs))

    coeffs = OrderedDict([(comp, c) for c, comp in zz])
    coeffs[products[0]] = n
    if verbosity >= 2:
        print("Debug\nCoeffs:{}".format(coeffs))
    show = []
    if verbosity >= 0:
        for compos in [reactants, products]:
            sh = ["{}{}".format(coeffs[comp], comp.iupac_formula) for comp in compos]
            sh = " + ".join(sh)
            show.append(sh)

        show = " ---> ".join(show)
        print(show)
    reac_coeffs = OrderedDict([(c, coeffs[c]) for c in reactants])
    prod_coeffs = OrderedDict([(c, coeffs[c]) for c in products])

    reac_sims_coeffs = OrderedDict([(org_r, coeffs[r]) for org_r, r in zip(org_reactants, reactants)])
    prod_sims_coeffs = OrderedDict([(org_p, coeffs[p]) for org_p, p in zip(org_products, products)])
    return (reac_coeffs, prod_coeffs), (reac_sims_coeffs, prod_sims_coeffs)


def equation_balancer_v1(reactants: list,
                      products: list,
                      verbosity: int = 1,
                      decimal_tolerance: int = 4,
                      depence_check: bool = False,
                      tol_zero: float = None ):
    '''Balances a given chemical equation,
    input parameters:
        reactants, products, verbosity(default=1)'''

    def fill_inmatrix(matrix, compounds: list, elsindex: np.asarray):
        for col, reacs in enumerate(compounds):
            for els in reacs:
                row = (elsindex == els.name).nonzero()
                assert len(row) == 1 and len(row[0]) == 1
                row = row[0][0]
                quan = reacs[els]
                matrix[row][col] = quan
        if verbosity >= 2:
            print("Debug:mode\n{Matrix}")
            print(matrix)
        return matrix

    org_reactants = reactants
    org_products = products

    reactants = [Pycomp(r) if isinstance(r, str) else r for r in reactants]
    products = [Pycomp(p) if isinstance(p, str) else p for p in products]

    lhs_els = np.unique([el.name for r in reactants for el in r])
    lhs_els.sort()
    rhs_els = np.unique([el.name for p in products for el in p])
    rhs_els.sort()

    if verbosity >= 2:
        print("debug:\nreactant_els:{}\nproduct_els:{}".format(lhs_els,rhs_els))

    if not np.all(lhs_els == rhs_els):
        raise ValueError("Reactants and products do not contain equal type of number of unique  elements")

    n = 1
    a = np.zeros((len(lhs_els), len(reactants)))
    b = np.zeros((len(rhs_els), len(products)))
    # fill_in
    a = fill_inmatrix(matrix=a, compounds=reactants, elsindex=lhs_els)
    b = fill_inmatrix(matrix=b, compounds=products, elsindex=rhs_els)
    orgb = b.copy()
    orga = a.copy()

    ab = -1 * b[:, 1:]

    try:
        a = np.hstack((a, ab))
    except ValueError:
        a = np.hstack((a, ab[np.newaxis].T))

    b = b[:, 0]
    assert b is not None

    if verbosity >= 2:
        print("Debug:\nSolving the MAtrix:\n{}X={}".format(a, b))

    if depence_check:
        n_org = len(a)
        rank = np.linalg.matrix_rank(a, tol=tol_zero)
        if a.shape[0] == a.shape[-1]:
            if n_org != rank:
                raise LinearlydependentMatrix("Matrix is not linearly-independent:{}".format(a))
        else:
            pass

    # Now Solve
    try:
        coeffs = solve(a, b)
    except ValueError:
        rows, columns = a.shape
        if rows >= columns:
            diff = abs(columns - rows)
            if depence_check:
                if verbosity >= 2:
                    print("Debug:\ndiff:{}".format(diff))
                if rank == rows or rank >= abs(rows-diff):
                    pass
                else:
                    raise LinearlydependentMatrix("Matrix is not linearly-independent:{}".format(a))

            combined_equation = np.hstack((a, b[:, np.newaxis]))

            if verbosity >= 2:
                print("Debug mode:\nFull equation {}".format(combined_equation))

            allsolutions = []
            for combin in itcomb(combined_equation, diff + 1):
                combin = np.asarray(combin)
                if verbosity >= 2:
                    print("Debug:\n{}\na={}\nb={}".format(combin, combin[:, :-1], combin[:, -1]))
                try:
                    coeffs = solve(combin[:, :-1], combin[:, -1])
                    allsolutions.append(coeffs)
                except LinAlgError:
                    continue
            if verbosity >= 2:
                print("Debug:\n solutions:{}".format(allsolutions))
            if len(allsolutions) < diff + 1:
                raise RuntimeError("Can't balance the equation")
            elif not all([np.allclose(allsolutions[0], b) for b in allsolutions]):
                raise RuntimeError("Can't balance the equation")

            coeffs = allsolutions[0]

        #     coeffs = solve(a[:diff], b[:diff])
        #     coeffs2 = solve(a[-diff:], b[-diff:])
        #     if verbosity >= 2:
        #         print("Debug:\n{}".format(coeffs))
        #         print("Coeffs_2:{}".format(coeffs2))
        #     if not np.allclose(coeffs, coeffs2):
        #         raise RuntimeError("Can't balance the equation")
        else:
            raise NotImplementedError("Unable to handle the matrix operations")

    except LinAlgError:  # non-singular Matrix
        coeffs = lstsq(a, b)[0]
        print(coeffs)

    #            try:
    #                coeffs,_,_,_ = lstsq(a,b)
    #                coeffs = np.around(coeffs, 3)
    #           except LinAlgError:
    #               raise LinAlgError("Couldnt balance the equation")
    try:
        factor = np.gcd.reduce(coeffs)
    except Exception:
        factor = 1

    coeffs = coeffs / factor
    n = n / factor
    if verbosity >= 2:
        print("Debug:\n{}".format(np.absolute(coeffs)))

    coeffs = np.around(coeffs, decimal_tolerance)

    if (coeffs == 0).any(-1):
        raise RuntimeError("One of the coefficient is zero: {}".format(coeffs))
    zz = zip(coeffs, reactants + products[1:])
    if verbosity:
        print("Debug:\n{}".format(coeffs))

    coeffs = OrderedDict([(comp, c) for c, comp in zz])
    coeffs[products[0]] = n
    if verbosity >= 2:
        print("Debug\nCoeffs:{}".format(coeffs))
    show = []
    if verbosity >= 0:
        for compos in [reactants, products]:
            sh = ["{}{}".format(coeffs[comp], comp.iupac_formula) for comp in compos]
            sh = " + ".join(sh)
            show.append(sh)

        show = " ---> ".join(show)
        print(show)
    reac_coeffs = OrderedDict([(c, coeffs[c]) for c in reactants])
    prod_coeffs = OrderedDict([(c, coeffs[c]) for c in products])

    reac_sims_coeffs = OrderedDict([(org_r, coeffs[r]) for org_r, r in zip(org_reactants, reactants)])
    prod_sims_coeffs = OrderedDict([(org_p, coeffs[p]) for org_p, p in zip(org_products, products)])
    return (reac_coeffs, prod_coeffs), (reac_sims_coeffs, prod_sims_coeffs)