from ase.spacegroup.symmetrize import FixSymmetry as FixSymmetryase


class FixSymmetry(FixSymmetryase):

    def todict(self):
        dcrt = {"name": FixSymmetryase.__name__}
        dcrt.update(self.__dict__)
        return dcrt

    def to_dict(self):
        return self.todict()

    def as_dict(self):
        dct = {"@module": self.__class__.__module__, "@class": self.__class__.__name__}
        dct.update(self.todict())
        return dct

    @classmethod
    def from_dict(cls, dct):
        from monty.json import MontyDecoder
        dct = MontyDecoder().process_decoded(dct)
        obj = cls(**dct)
        return obj