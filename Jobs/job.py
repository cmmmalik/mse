from ase.atoms import Atoms
from ase.io import read, write
from ase.db import connect as asedbconnect
from ase.db.core import Database as dbCore
import os
import sys
from pathlib import Path as pathlibPath
import subprocess
import socket
import traceback
from typing import Any
import warnings


from mse.system import directory
from Database.ASE_io import Gpawjobdb_read
import HPCtools.hpc_tools3 as HPC

# TODO: add function for cleaning the contents of the directory instead of the directory itself.


class Corejob:

    _keywords = {"init": "initialized",
                 "fin": "finalized",
                 "run": "running",
                 "rm": "removed",
                 "res": "reset"}

    def __init__(self, name, working_directory=None):
        self._name = name

        if not working_directory:
            working_directory = os.getcwd()

        self.working_directory = working_directory
        self._backup_folder = None
        self._status = None
        self._host_directory = None
        self._host_address = None
        self._local_name = None
        self._autofetch = None
        self._id = None

    def __repr__(self):
        st = "working_directory:{wrk}".format(wrk=self.working_directory)
        if self.status:
            st = st + "status:{st}".format(st=self.status)

        return "{0}({1})".format(self.__class__.__name__, st)

    def __str__(self):
        return self.__repr__()

    def __getstate__(self): # needed to pickle
        dct = self.__dict__.copy()
        dct["_keywords"] = self._keywords # __dict__ does not class (inherent) attributes like _keywords
        #TODO: do this in a better way: create a separate attribute of _keywords
        return dct

    def __setstate__(self, state):
        self.__dict__.update(state)

    @property
    def autofetch(self):
        return self._autofetch

    @autofetch.setter
    def autofetch(self, value:bool):
        if isinstance(value, bool):
            self._autofetch = value
        else:
            raise ValueError(f"Expected an instance of {bool} , but got {type(value)}")

    @property
    def name(self):
        return self._name

    @property
    def host_directory(self):
        return self._host_directory

    @property
    def local_name(self):
        return self._local_name

    @property
    def host_address(self):
        return self._host_address

    @property
    def Id(self):
        return self._id

    @property
    def backup_folder(self):
        return self._backup_folder

    def initialize(self):
        self.create_dir()

    def _setinitialized(self):
        self._setstatus(self._keywords["init"])

    def _setfinished(self):
        self._setstatus(self._keywords["fin"])

    def _setrunning(self):
        self._setstatus(self._keywords["run"])

    def _setremoved(self):
        self._setstatus(self._keywords["rm"])

    def _setreset(self):
        self._setstatus(self._keywords["res"])

    def create_dir(self):
        try:
            os.mkdir(self.working_directory)
        except FileNotFoundError:
            os.makedirs(self.working_directory)

    def backup_dir(self):
        if not os.path.exists(self.working_directory):
            warnings.warn("THe folder does not exist, so leaving it untouched")
            return
        self._backup_folder = self.get_backup_name()
        os.mkdir(self.backup_folder)

    @staticmethod
    def _backup_name():
        return ".Backup"

    def get_backup_name(self):
        c = 0
        backup = self._backup_name()
        while True:
            fullbackup = os.path.join(self.working_directory, backup)
            if not os.path.exists(fullbackup):
                return fullbackup
            c += 1
            backup = "{}_{}".format(backup, c)

    @property
    def keywords(self):
        return self._keywords

    def remove(self):
        directory.Directory.remove_folder(self.working_directory)
        self._setremoved()

    def remove_force(self, safe=True):
        if safe:
            self.reset()
        else:
            warnings.warn(f"The {self.working_directory} will be permanently deleted whether it's empty or not")
            ui = input("Enter yes(no):")
            if ui in ["yes", "y"]:
                directory.Directory.delete(self.working_directory)
        self._setremoved()

    def reset(self):
        #clean directory
        status = directory.Directory.clean_folder(self.working_directory)
        if status:
            self._setreset()
        else:
            warnings.warn("Couldn't reset the working directory:{}".format(self.working_directory))
            print("Manually remove the contents")

    @property
    def status(self):
        return self._status

    def _setstatus(self, value: str):
        self._status = value

    def run_check(self):
        if self.status != self._keywords["init"]:
            raise RuntimeError(f"The job is not prepared/initialised properly, Found status: {self.status}")

    @staticmethod
    def read_path(path: str, name=None, ):
        if os.path.exists(path):
            return Corejob(name=name, working_directory=path)
        else:
            raise ValueError(f"Path '{path}' does not exist")

    @property
    def working_directory(self):
        return self._working_directory

    @working_directory.setter
    def working_directory(self, peth):
        self._working_directory = pathlibPath(peth).expanduser().resolve()


class ASEjob(Corejob):

    def __init__(self, name, working_directory=None, atoms=None):
        super(ASEjob, self).__init__(name, working_directory=working_directory)

        self.atoms = atoms
        self._hpc = None

    def __repr__(self):

        st = super(ASEjob, self).__repr__()
        st = st.strip(")")

        if self.atoms:
            st += ",atoms:{atoms}".format(atoms=self.atoms)

        st = st + ")"
        return st

    def __str__(self):
        return self.__repr__()

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, value: Atoms):
        if isinstance(value, Atoms):
            self._atoms = value
        elif not value:
            self._atoms = None
        else:
            raise ValueError(f"Expected an instance of {Atoms}, recieved {type(value)}")

    @staticmethod
    def atoms_frompath(filename,
                       index: Any = None,
                       format: str = None,
                       **kwargs):

        atoms = read(filename, index=index, format=format, **kwargs)
        return atoms

    def structure_write(self, filename, **kwargs):
        write(filename, self.atoms, **kwargs)

    @classmethod
    def read_path(cls, name=None, working_directory=None):
        obj = cls.__init__(name=name, working_directory=working_directory)
        obj.atoms = obj.atoms_frompath()
        return obj

    @staticmethod
    def ase_row_db(db, index: int = -1, verbosity: int = 0):
        rows = Gpawjobdb_read(db=db).R
        if verbosity >= 1:
            print(rows)
            print("Found rows:")
            for r in rows:
                print("{}{}".format(r.id, r))

        if len(rows) > 1:
            warnings.warn("Found more than one rows in the database {}"
                          "\nSelecting the last index row\nUse 'index' to change the read row".format(db))
        rows = rows[index]
        return rows

    @staticmethod
    def _smartwriteasedb(db: str or dbCore, row, data: dict = {}, keys: dict = {}, **kwargs):

        if isinstance(db, dbCore):
            db.write(row, data=data, key_value_pairs=keys, **kwargs)
        else:
            with asedbconnect(db) as mydb:
                mydb.write(row, data=data, **keys)

    @property
    def hpc(self):
        return self._hpc

    @hpc.setter
    def hpc(self, value):
        if isinstance(value, HPCjob):
            self._hpc = value
        else:
            raise TypeError(f"Expected an instance of {HPCjob}, got {type(value)}")

        try: # also set the local machine address
            self._local_name = self.hpc.grab_localname()
        except Exception as e:
            print(f"***Encountered Exception {e}")
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)

    def initialize_hpc(self,
                       server: str = "lcluster13",
                       tipe: str = None,
                       verbosity: int = 0):
        hpc = HPCjob(jobname=self.name,
                     working_directory=self.working_directory,
                     server=server,
                     tipe=tipe,
                     verbosity=verbosity)
        self.hpc = hpc
        return hpc


class HPCjob:

    def __init__(self,
                 jobname: str,
                 working_directory: str = None,
                 server: str = "lcluster",
                 tipe: str = None,
                 verbosity: int = 0):

        self._init_pre(jobname=jobname, working_directory=working_directory)

        self._server = HPC.HPCMain(server=server,
                                   tipe=tipe,
                                   jobname=getattr(self, "jobname", "HPCtoolsrun_v2"),
                                   localdir=self.working_directory,
                                   verbosity=verbosity)

    def _init_pre(self,
                  jobname: str,
                  working_directory: str=None):

        self._jobname = jobname

        if working_directory:
            self.working_directory = working_directory
        else:
            self._working_directory = os.getcwd()

        self._main_directory = os.getcwd()

    def __del__(self):
        if os.getcwd() != self.main_directory:
            try:
                os.chdir(self.main_directory)
            except:
                pass

    def __repr__(self):
        st = ""
        if self.jobname:
            st += "jobname:{},".format(self.jobname)
        if self.working_directory:
            st += "working_directory:{},".format(self.working_directory)
        if self.main_directory:
            st += "main_directory:{},".format(self.main_directory)
        if self.server:
            st += "server:{}".format(self.server)

        return "{0}({1})".format(self.__class__.__name__, st)

    def __str__(self):
        return self.__repr__()

    @property
    def jobname(self):
        return self._jobname

    @property
    def working_directory(self):
        return self._working_directory

    @working_directory.setter
    def working_directory(self, peth):
        self._working_directory = pathlibPath(peth).expanduser().resolve()

    @property
    def main_directory(self):
        return self._main_directory

    @property
    def server(self):
        return self._server

    def prepare(self,
                random_folder = True,
                folder_name: bool = False,
                **kwargs):
        self.server.prepare(random_folder=random_folder, folder_name=folder_name, **kwargs)

    def submit(self):
        try:
            self.server.submit()
        except Exception:
            raise Exception("Couldn't submit the {}".format(self.jobname))
        finally:
            # revert back to main directory
            os.chdir(self.main_directory)

    def run_hpc(self, code: str = None, **kwargs):
        self.prepare(code=code, **kwargs)
        self.submit()

    def fetch_hpc(self, **kwargs):
        try:
            self.server.fetch(**kwargs)
        except Exception:
            raise Exception("Couldn't fetch the job, see traceback")
        finally:
            os.chdir(self.main_directory)

    @staticmethod
    def grab_localfullname():
        proc = subprocess.run(["hostname", "-A"], check=True, capture_output=True, text=True)
        name = proc.stdout
        # grab the top one
        name = name[-1]
        return name

    @staticmethod
    def grab_localname():
        return socket.gethostname()

    @classmethod
    def read_hpc(cls,
                 jobname: str = None,
                 working_directory: str = None,
                 iid: str or int = None,
                 verbose: int = 0,
                 jobid: str or int = None,
                 **kwargs):

        # jobname is not important while reading back
        cls = cls.__new__(cls) # conerting class to instance
        cls._init_pre(jobname=jobname, working_directory=working_directory)
        cls._server = HPC.HPCMain.read_hpc(iid=iid, verbosity=verbose, jobid=jobid, working_directory=working_directory, **kwargs)
        return cls

    @classmethod
    def from_server_hpc(cls,
                        hostname: str,
                        jobname: str = None,
                        working_directory: str = None,
                        iid: str or int = None,
                        jobid: str or int = None,
                        hpcdir: str = None,
                        localdir: str = None,
                        verbosity: int = 0):

        cls = cls.__new__(cls) # convert class to instance
        cls._init_pre(jobname=jobname, working_directory=working_directory)
        cls._server = HPC.HPCMain.recreate_basic_hpc(hostname=hostname,
                                                     jobname=jobname,
                                                     jobid=jobid,
                                                     iid=iid,
                                                     workdir=hpcdir,
                                                     localdir=localdir,
                                                     verbosity=verbosity)
        return cls