from ase.atoms import Atoms

from mse.workflows.base import Baseworkflow
from mse.optimizer.optimizer import Moptimize


class Workflow_relaxation(Baseworkflow):

    def __init__(self,
                 atoms: Atoms,
                 working_directory: str,
                 reltype: list or str = ["ions", "cell", "full"],
                 verbosity: int = 1,
                 dircheck: bool = True):

        super(Workflow_relaxation, self).__init__(atoms=atoms,
                                                  working_directory=working_directory,
                                                  calculator_type="gpaw",
                                                  verbosity=verbosity,
                                                  dircheck=dircheck)
        self._runtype = "relax"
        if not isinstance(reltype, (list or tuple)):
            reltype = [reltype]
        self._reltype = reltype


    @property
    def reltype(self):
        return self._reltype

    @property
    def runtype(self):
        return self._runtype

    def set_optimizer(self):
        self.job.newrunscheme = Moptimize
        self.job.relax_inputs["reltype"] = self.reltype

    def initialize_job(self,
                       name: str,
                       modeinps: dict = None,
                       calcinps: dict = None,
                       **kwargs):

        kwargs.pop("run_type", None)
#        reltype = kwargs.get("relaxargs", None)
#        try:
#           self._reltype = reltype.pop("reltype", self.reltype)
#        except AttributeError:
#            pass
        super(Workflow_relaxation, self).initialize_job(name=name,
                                                        run_type="relax",
                                                        modeinps=modeinps,
                                                        calcinps=calcinps,
                                                        **kwargs)
        self.set_optimizer()

        return self.job
