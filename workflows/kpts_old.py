import numpy as np
import os
import warnings

from ..system.directory import SpDirectory

#workflow for performing k-mesh convergence--old


class Workflow_kp:

    def __init__(self,
                 el,
                 poscarpath,
                 extrasfolder_tags=None,
                 encut=1200,
                 kplimits=[4, 9],
                 outfile="run.py",
                 safe=True):

        #TODo: all inputs should be safed in a self.input dictionary

        self.el = el
        self.poscarpath = poscarpath
        self.encut = encut
        self.kplimits = kplimits
        self.outfile = outfile

        if extrasfolder_tags:
            self.folder = el + extrasfolder_tags
        else:
            self.folder = el

        self.ins = self._initinspect(safe)

        if not self.ins:
            print("I will not proceed")

        self.subdirs = None
        self.kp = None

    def _initinspect(self, safe):

        exists = os.path.exists(self.folder)

        if not exists:
            return True
        if not SpDirectory.folder_empty(self.folder):
            warnings.warn("The folder is not empty")
            if safe:
                return False
            else:
                print("Even though, you have opted to proceed")

        return True

    def autoflow(self):

        if not self.ins:
            print("I will not proceed")
            return

        # build up subdirectores

        self.initialize()

        print(f"Creating directories '{self.subdirs}'")

        self.create_subdirs()
        self.setup_gpaw()

        # Show_Tree structures
        SpDirectory.display_tree_structure(self.folder)

    def initialize(self):

        kp = np.arange(self.kplimits[0], self.kplimits[1])
        self.kp = kp
        subnames = [os.path.join("%s" % self.folder, "%s.%s" % (self.encut, k)) for k in kp]
        self.subdirs = subnames

        print(f"The subdirs will be created are: {self.subdirs}")

    def create_subdirs(self):

        SpDirectory.buildup_kp(self.subdirs)

    def setup_gpaw(self, poscarname="POSCAR"):

        SpDirectory.copy_poscars(self.subdirs, self.poscarpath)
        SpDirectory.write_gpawscript(self.subdirs, encut=self.encut, outfile=self.outfile, poscarname=poscarname)
        SpDirectory.display_tree_structure(self.folder)

    @staticmethod
    def hpc_submit(folders,
                   hpcname="lcluster",
                   walltime="3:00:00",
                   code="gpaw"):

        # ToDo: put HPC in a separate class,

        import hpc_tools3 as HPC

        nosuccess = []
        server = HPC.HPCtools(hpcname)
        # have to jump into directory#
        maindirec = os.getcwd()

        for f in folders:
            os.chdir(f)
            server.params["walltime"] = walltime
            server.jobname = "%s" % f
            server.prepare_job(code=code)
            try:

                server.submit_job()

            except Exception:
                print("I couldn't submit-something is wrong")
                nosuccess.append(f)

            os.chdir(maindirec)
        if nosuccess:
            return nosuccess

        return True

    @staticmethod
    def hpc_fetch(folders,
                  hpcname="lcluster"):

        import hpc_tools3 as HPC
        server = HPC.HPCtools(hpcname)

        nosucess = []
        maindirec = os.getcwd()

        for f in folders:
            os.chdir(f)
            try:

                server.fetch_job()
            except Exception:
                print("I couldn't fetch the job")
                nosucess.append(f)
            os.chdir(maindirec)

        if nosucess:
            return nosucess

        return True

# should be part of workflow, a naive way of finding the value under which kp is converged


def get_conv_index(energies=None, diff_energies=None, tolerance=0.5):
    assert energies or diff_energies
    if energies:
        delta = [abs(i - energies[-1]) * 1000 for i in energies]
    else:
        print("Make sure the units of tolerance and values in the arrays 'energies' or 'diff_energies' are same")
        delta = diff_energies

    index = np.argmin(delta[:-1])

    if delta[index] < tolerance:
        return index

    print("Tolerance is too high or the energies are not converged")
    return None

