from pymatgen.ext.matproj import MPRester


def get_entries(formula,
                mpkey,
                compatible_only=False,
                sort_by_e_above_hull=True,
                inc_structure=None,
                **kwargs):

    with MPRester(mpkey) as mp:
        entries = mp.get_entries(formula,
                                compatible_only=compatible_only,
                                inc_structure=inc_structure,
                                sort_by_e_above_hull=sort_by_e_above_hull,
                                **kwargs)
    return entries


class SmartMPRester:

    def __init__(self, mpkey):
        self._mpkey = mpkey
        self._connection = None

    @property
    def mpkey(self):
        return self._mpkey

    @property
    def connection(self):
        return self._connection

    def __del__(self):
        if self.connection:
            try:
                self.close()
            except Exception:
                pass

    def connect(self):
        self._connection = MPRester(api_key=self._mpkey)

    def close(self):
        self.connection.session.close()

    def get_entries(self, *args, **kwargs):
        return self.connection.get_entries(*args, **kwargs)