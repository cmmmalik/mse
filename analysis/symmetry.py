
def get_symmetry_kp_gamma(symmetry, verbosity: int = 1):
    number = symmetry["number"]
    number = int(number)
    if verbosity >= 1:
        print("Spacegroup number: {}".format(number))

    if number >= 168:  # then we have hexagonal spacegroups
        return {"gamma": True}
    else:
        return dict()