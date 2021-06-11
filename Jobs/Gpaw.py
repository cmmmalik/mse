from ase.db import connect, core as dbCore
import os
import json
import subprocess
import sys
import traceback
import warnings

from Database.ase_db import Savedb
from mse.Jobs.job import ASEjob, HPCjob
from mse.io.scripts import InterfaceGPAW
import mse.system.directory as Simdirec
from HPCtools_v2.hpc_tools3 import directorychange

# TODO: Implement JOBsets
# TODO: Implement job copy method
# TODO: Implement is_converged method, if calc is present then it is converged in case of static calculation,
#   otherwise,check for the forces
# TODO: Child jobs similar to pyiron
# TODO: Move HPC to the Core job or subclass Gpaw with HPCjob. Depending upon the design of the project.
# TODO: Dump more information about job into the inputs.json file
# TODO: Shift some methods to Core job. or ASEJobs
# TODO: Add a method to re_run the job
# TODO: Add .copy method
# TODO: remove outputs upon cleaning/restarting/reruning the job
# TODO: separate the attributes needed for using InterfaceGPAW writer.


class Gpaw(ASEjob):

    _keywords = {"inpjson": "inputs.json", "pyname": "run.py"}

    def __init__(self, name, working_directory=None, atoms=None):

        super(Gpaw, self).__init__(name=name, working_directory=working_directory, atoms=atoms)

        self.inputs = dict()
        self._defaults_inputs()

        self._outputs = {}
        self._gpaw_reader = None

        self._pyname = self._keywords["pyname"]
        self._hpc = None

        self._keywords.update(super(ASEjob, self)._keywords)

        self._caltype = None
        self._mtype = None
        self._code = "gpaw"

  #  def __getattr__(self, item):
   #     if self._outputs:
    #        try:
   #             return self._outputs[item]
   #         except KeyError:
   #             raise AttributeError(f"instance {self.__class__.__name__} does not have attribute {item}")

    def __repr__(self): # TODO update this with more relevant attributes
        return super(Gpaw, self).__repr__()

    def __str__(self):
        return self.__repr__()

    def _defaults_inputs(self):
        self.inputs = {"encut": 500, "kden": 4, "xc": "PBE", "poscarname": "POSCAR", "outputtxt": None}

    @directorychange
    def getoutputs(self):

        gpaw_reader = Savedb.gpawcollect(txtfile=self.inputs["outputtxt"],
                                         pyfile=self._keywords["pyname"],
                                         poscarfile=self.inputs.get("poscarname"),
                                         )
        if self._outputs:
            raise RuntimeError(f"A possible bug, attribute outputs: {self.outputs} is not empty")

        self._gpaw_reader = gpaw_reader
        self._outputs["atoms"] = gpaw_reader.atoms
        self._outputs["init_atoms"] = gpaw_reader.init_atoms
        self._outputs["inps"] = gpaw_reader.calc
        self.atoms = gpaw_reader.atoms

    def _static(self, mtype="bulk"): # TODO: implement other types of relaxations.
        if mtype == "bulk":
                self.inputs["outputtxt"] = "static.txt"
                writer = InterfaceGPAW() # TODO: very basic, need to extend and develop this interface for writing script files
                writer.static_script_bulk(outfile=os.path.join(self.working_directory, self.pyname), **self.inputs)

    @property
    def caltype(self):
        return self._caltype

    @caltype.setter
    def caltype(self, value:str):
        if not self.status:
            self._caltype = value
        else:
            raise RuntimeError(f"The job has been already initialized , can't set the attribute 'caltype' to '{value}'")

    def fetch_hpc(self):
        self.hpc.fetch_hpc()
        self._setfinished()

    @property
    def gpaw_reader(self):
        return self._gpaw_reader

    def grabhpc(self, server="lcluster"):
        return HPCjob(jobname=self.name, working_directory=self.working_directory, server=server)

    @property
    def hpc(self):
        return self._hpc

    @hpc.setter
    def hpc(self, value: HPCjob):
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

    def initialize(self):
        super(ASEjob, self).initialize()
        super(Gpaw, self).structure_write(
            filename=os.path.join(self.working_directory, getattr(self.inputs,"poscarname", "POSCAR")),
            format="vasp",
            vasp5=True)

    def is_notconverged(self): # Use this if you want to know job is not converged
        return False if self.atoms.calc else True

    @property
    def keywords(self):
        return self._keywords

    @property
    def mtype(self):
        return self._mtype

    @mtype.setter
    def mtype(self, value):
        if not self.status:
            self._mtype = value
        else:
            raise RuntimeError(f"The job has been already initialized , can't set the attribute 'caltype' to '{value}'")

    @property
    def outputs(self):
        return self._outputs

    def prepare(self, caltype="static", mtype="bulk"):
        self.caltype = caltype
        self.mtype = mtype

        if mtype != "bulk":
            raise NotImplementedError(f"Currently only supports 'mtype':'bulk' , not {mtype}")

        if caltype == "static":
            self._static(mtype=mtype)

        self.saveinputs(filename=os.path.join(self.working_directory, self._keywords["inpjson"]))
        self._setinitialized()

    @property
    def pyname(self):
        return self._pyname

    @classmethod
    def read_path(cls,
                  path: str,
                  name: str = None,
                  outputtxt: str = None):

        if not os.path.exists(path):
            raise ValueError(f"Path '{path}' does not exist")

        json_loc = os.path.join(path, cls._keywords["inpjson"])
        if os.path.exists(json_loc):
            with open(json_loc, "r") as jsonf:
                inputs = json.load(jsonf)

        else:
            warnings.warn("At the moment Json file must be present,"
                          "'.inputs' will be set to None"
                          "Only use '.outputs' attribute", DeprecationWarning)

            inputs = {}

        obj = Gpaw(name=name, working_directory=path)

        if inputs:
            obj.set_inputs(**inputs)
            obj._mtype = obj.inputs["outputtxt"]
            # TODO: Work around for setting job.mtype, better would be to save mtype also in inputs,json
        else:
            obj.inputs = inputs
            if not outputtxt:
                outputtxt = Simdirec.Directory.findglob_file("*.txt", directory=path)
                assert len(outputtxt) == 1
                outputtxt = outputtxt[0]

        if outputtxt:
            obj.inputs["outputtxt"] = outputtxt

        obj.getoutputs()
        return obj

    def run(self, wait=False):
        if self.hpc:
            self.submit_hpc()
            self._setrunning()
        else:
            print("Running locally, HPC server not found")
            if not wait:
                self.run_local()
            else:
                self.run_local_wait()

    def run_local(self):
        super(ASEjob, self).run_check()
        p = subprocess.Popen(["gpaw", "python", "%s" %os.path.join(self.working_directory, self.pyname)],
                           stdout=None)
        self._setrunning()
        return p

    def run_local_wait(self):
        super(ASEjob, self).run_check()
        p = self.run_local()
        p.wait()
        self._setfinished()

    def saveinputs(self, filename="inputs.json"):
        with open(filename, "w") as f: # it will overwrite the already existing file
            json.dump(self.inputs, f)

    def save_todatabase(self, db: str or dbCore):
        self._gpaw_reader.get_keys_asedb()
        data = self._gpaw_reader.output_data_asedb()
        keys = self._gpaw_reader.keys_parsedb()

        try:
            # if not db.closed:
            db.write(self.outputs["atoms"], data=data, **keys)
        except (AttributeError, IOError) as e:
            with connect(db) as mydb:
                mydb.write(self.outputs["atoms"], data=data, **keys)
        except Exception:
            raise Exception(f"Couldn't write into the database: {db}, Follow the traceback")

    def set_inputs(self, **kwargs):
        self.inputs.update(kwargs)

    def submit_hpc(self):
        super(ASEjob, self).run_check()
        self.hpc.run_hpc(code="gpaw")
        self._setrunning()

    def detach_hpc(self):
        hpc = self.hpc
        self._hpc = None
        return hpc