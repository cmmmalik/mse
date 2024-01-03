import os
from pathlib import Path
from .job import HPCjob


#...At the moment as simple wrapper....

class Hpcarray:

    def __init__(self,
                 jobname:str,
                 working_directory:str=None,
                 ):

        self.name = jobname
        self.working_directory = working_directory if working_directory else os.getcwd()
        self._jobs = None
        self._hpc = None
        self._type = None

    @property
    def jobs(self):
        return self._jobs

    @property
    def type(self):
        return self._type

    @jobs.setter
    def jobs(self, js:list or tuple):
        assert [Path(j.working_directory).parent == Path(self.working_directory).name for j in js]
        self._type = getattr(js[0], "calc")
        self._jobs = js

    def set_jobs(self, jobs,):
        self.jobs = jobs

    def set_type(self, type):
        self._type = type

    def initiate(self,
                 server: str = None,
                 tipe: str = None):

        from HPCtools.hpc_arrays import HPCarray
        self._server = HPCarray(server=server,
                                tipe=tipe,
                                jobname=self.name,
                                verbosity=0,
                                localdir=self.working_directory,
                                )

        self._server._initial_check()
        self._server.write_submit(code=self.type.lower(),)





