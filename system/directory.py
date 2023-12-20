from colorama import Style, Fore, init
from functools import partial
import fnmatch
import glob
import os
import shutil
import subprocess
import sys
import warnings


init(autoreset=True)


def sub_runspecific():

    return partial(subprocess.run, check=True, capture_output=True, text=True)


def subrun(cmd, **kwargs):
    pfunction = sub_runspecific()
    return pfunction(cmd, **kwargs)


class Directory:
    @classmethod
    def mkdir(cls, path):

        cmd = "mkdir -p %s" %path
        cmd = cmd.split()
        proc = cls.runcmd(cmd)
        return proc

    @classmethod
    def cp(cls, pfrom, pto):

        print("Copy the files %s to %s" %(pfrom, pto))
        cmd = "cp %s %s" %(pfrom, pto)
        cmd = cmd.split()
        proc = cls.runcmd(cmd)
        return proc

    @staticmethod
    def remove_folder(directory):
        try:
            os.rmdir(directory)
        except FileNotFoundError:
            print(f"Working directory {directory} does not exist")
        except OSError:
            print(f"Working directory {directory} is not empty")
            print(f"Either remove the folder manually or use 'delete' with safe=False")

    @classmethod
    def clean_folder(cls, directory, force=False):
        contents = os.listdir(directory)
        for i in contents:
            try:
                os.remove(os.path.join(directory, i))
            except OSError:
                if force:
                    cls.delete(os.path.join(directory, i))
                warnings.warn("Can not remove the folder {}".format(i))
        else:
            return True

        return False

    @staticmethod
    def delete(directory):
        try:
            shutil.rmtree(directory)
        except Exception as e:
            print("The directory couldn't be removed")
            print(e)
            raise e

    @classmethod
    def scp(cls, pfrom, pto):

        cmd = "scp %s %s" % (pfrom, pto)
        cmd = cmd.split()
        proc = cls.runcmd(cmd)
        return proc

    @classmethod
    def mv(cls, pfrom, pto):

        cmd = "mv %s %s" %(pfrom, pto)
        cmd = cmd.split()
        proc = cls.run(cmd)
        return proc

    @classmethod
    def rsync(cls, pfrom, pto, flags="-uavzh", dry=False):

        command = "rsync {} {} {} ".format(flags, pfrom, pto)
        if dry:
            command += " -n"

        command = command.split()
        proc = cls.runcmd(command)
        return proc

    @staticmethod
    def display_tree_structure(path="."):

        """Function to display the subfolders as tree structure """

        for root, dirs, files in os.walk(path):

            level = root.replace(path, "").count(os.sep)
            indent = ' ' *4* (level)
            print(Fore.YELLOW + Style.BRIGHT + '{}{}/'.format(indent, os.path.basename(root)))

            subindent = ' ' *4*(level + 1)

            for f in files:
                print(Style.BRIGHT+Fore.RED+ '{}{}'.format(subindent, f))

    @staticmethod
    def runcmd(cmd, **kwargs):
        proc = subrun(cmd, **kwargs)
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)

        return proc


    @staticmethod
    def findglob_file(query, directory):

        """

        @param query: (str) keyword to match with a file
        @param directory: (str) path in which th query(-containing) file will be searched
        """
        if not isinstance(query, str):
            raise ValueError(f"query '{query}' must be an instance of a string '{str}'")

        return glob.glob(os.path.join(directory, query)) # ONly non-hidden files

    @staticmethod
    def findfnmatch_file(directory:str, allquery=None, anyquery=None,case_sensitive=False):

        assert allquery or anyquery

        if allquery:
            cond = all
            query = allquery

        elif anyquery:
            cond = any
            query = anyquery

        contents = os.listdir(directory)

        if not contents:
            return None

        if isinstance(query, str):
            query = [query]

        Fnmatch = fnmatch.fnmatch

        if case_sensitive:
            Fnmatch = fnmatch.fnmatchcase

        return [os.path.join(directory,i) for i in contents if cond([Fnmatch(i, q) for q in query])]


class SpDirectory(Directory):

    @staticmethod
    def folder_empty(path):
        # currently checks for only non-hidden files
        files = os.listdir(path)
        files = [f for f in files if not f.startswith(".")]
        if files:
            return False

        return True

    @staticmethod
    def copy_poscars(folders, poscar):
        # poscars copy -store them in input database also!
        for k in folders:
            shutil.copy(poscar, k)

    @staticmethod
    def buildup_kp(folders):
        for k in folders:
            os.makedirs(k)


    @staticmethod
    def write_gpawscript(folders, encut, outfile="run.py", poscarname="POSCAR", sep="."):
        # k-points density (kden) and encut, and poscarname="POCSAR"
        from mse.io.scripts import InterfaceGPAW as itfgpaw
        writer_gpaw = itfgpaw(poscarname)

        for k in folders:
            kden = k.rsplit(sep, maxsplit=1)[-1]
            if not kden:
                raise RuntimeError(f" 'kden' is not present in the folders name")

            writer_gpaw.static_script_bulk(encut=encut, kden=kden, outfile=os.path.join(k, outfile))


class CD:

    """Manager for changing to the directories. Keep track of the old directories """

    def __init__(self):

        self.oldpaths = []
        self.newpath = None

    def enter(self, newpath):
        self.oldpaths.append(os.getcwd())
        os.chdir(newpath)
        self.newpath = newpath
        print("Jumped to directory {}".format(self.newpath))

    def exit(self):
        try:
            os.chdir(self.oldpaths[0])
        except IndexError as ex:
            pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.exit()

    def stepback(self):
        self.directjump(index=-1)

    def directjump(self, index):
        print(self.oldpaths)
        try:
            os.chdir(self.oldpaths[index])
            print("Back Jumped to directory {}".format(self.oldpaths[index]))
            self.oldpaths = self.oldpaths[:index]
        except IndexError:
            raise RuntimeError("Already at the root directory/ or index is incorrect")

    def __repr__(self):
        return "Current directory :" + os.getcwd()


def directorychange(func):

    def dectorator(ins, *args, **kwargs):

        cd = CD()
        cwd = os.getcwd()
        working_directory = kwargs.pop("working_directory", None)
        caution = kwargs.pop("caution", None)

        try:
            localdir = ins._localdir
        except AttributeError:
            localdir = getattr(ins, "_working_directory", working_directory)

        if not localdir:
            localdir = cwd
        print(localdir)
        if str(cwd) != str(localdir):
            try:
                cd.enter(localdir) # should give an error, if directory does not exist
            except OSError:
                print("Couldn't go into the directory, {} does not exist".format(localdir))
                if caution:
                    raise OSError
        try:
            out = func(ins, *args, **kwargs)
        except Exception as exc:
            raise exc

        finally:
            try:
                cd.exit()
            except IndexError:
                pass # means that we did not move at all
        return out

    return dectorator
