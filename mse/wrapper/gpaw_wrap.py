import pickle
import argparse
import warnings

from ase.parallel import paropen, parprint as aseparprint, world

from ase.db import connect as asedbconnect
from copy import deepcopy
from functools import partial
import traceback
import sys

# from mse.calculator.gpaw_calc import Gpaw

parprint = partial(aseparprint, flush=True) # flush the buffers
print = partial(print, flush=True)

# TODO: Add option of debugging, that will increase the verbosity


def generate_runcmd(inputfile,
                    save_calc: bool = True,
                    save_db: bool = False,
                    backup: bool = False,
                    remove: bool = True,
                    envcmd: bool = None):

    loc = "mse.wrapper.gpaw_wrap"
    runcmd = "mpirun gpaw python -m {}".format(loc)
    runcmd += " -pf {}".format(inputfile)
    runcmd += " --savecalc {}".format(save_calc)
    runcmd += " --savedb {}".format(save_db)
    runcmd += " --remove {0} --backup {1}".format(remove, backup)
    if envcmd is not None:
        runcmd += " -e {}".format(envcmd)
    return runcmd


def commandlineargs():

    parser = argparse.ArgumentParser()
    parser.add_argument("--unpicklef", "-pf",
                        help="File to be unpickled, which contains the job instance",
                        default="job.p",
                        type=str,
                        dest="filename")

    parser.add_argument("--savedb", "-s",
                        help="Save the atoms into database .db, can bt used in spite of saving the whole calculator",
                        dest="savedb",
                        default=False)

    parser.add_argument("--savecalc", "-C",
                        help="Save the calculator",
                        dest="savecalc",
                        default=True)

    parser.add_argument("--remove", "-r",
                        help="set True(or False) to remove(or not ) the folder, respectively",
                        default=False,
                        dest="remove")

    parser.add_argument("--backup", "-b",
                        help="Set to True(or False) to backup(or not) the local folder contents",
                        default=False,
                        dest="backup")

    parser.add_argument("--envcmd", "-e",
                        help="Extra environmental bash commands that shall be executed before invoking the fetcher e.g."
                             "'conda activate' etc. (Currently supports only a single command)",
                        type=str,
                        dest="envcmd")

    args = parser.parse_args()
    return args


def initialize_calc(job): #ToDo: this should be part of the job -----

    import gpaw

    CALC = getattr(gpaw, job.calc)
    MODE = getattr(gpaw, job.mode)
    inputs = deepcopy(job.inputs)

    calc = None
    modeargs = inputs["mode_args"]

    try:
        ecut = modeargs.pop("encut") # better to change this key to ecut
    except KeyError:
        warnings.warn("Encut was not provided", RuntimeWarning)
        ecut = None

    attach = inputs["calc_args"].pop("attach", None)

    if job.restart:
        from gpaw import restart
        atoms, calc = restart(filename="calc.gpw", **inputs["calc_args"])
        if job.atoms:
            warnings.warn("Atoms object was found in the job, it will be overwritten with atoms read from calc.gpw ")

        job.atoms = atoms
        # job.atoms.calc = calc(restart="calc.gpw", mode=mode(ecut=ecut, **modeargs), **inputs["calc_args"])
    if not calc:
        calc = CALC(mode=MODE(ecut=ecut, **modeargs), **inputs["calc_args"])

    if attach:
        calc.attach(calc.write, **attach)

    job.atoms.calc = calc


def read_jobinfo(jobinfo: str = "job.info"):

    target = None
    origin = None

    with paropen(jobinfo, "r") as f: # this will read on all the slaves
        for line in f:
            line = line.strip()
            if "ORIGIN" in line:
                origin = line.split(maxsplit=2)[1:]
            elif "TARGET" in line:
                target = line.split(maxsplit=2)[1:]
            elif "JOBID" in line:
                jobid = line.rsplit(maxsplit=1)[-1]
            elif "JOBNAME" in line:
                jobname = line.split(maxsplit=1)[-1]
        return origin, target, jobid, jobname


def autosendback(remove:bool, backup:bool, envcmd:str, verbosity:int=0):

# Host refers to non-local machine, present on the other end of ssh connection.
# Do not raise any exception, let the job pass through to completion.
    from mse.servertools.Slurm import putback_origin
    from HPCtools.hpc_tools3 import HPCMain

    origin, target, jobid, jobname = read_jobinfo(HPCMain.jobinfo)
    hostname, hostdir = origin
    locname, locdir = target
    parprint("Host (to which job is being fetched back):\n{0}, {1}".format(hostname, hostdir))

    if world.rank == 0: # use only master to send back

        try:
            back = putback_origin(jobname=jobname,
                                  host=hostname,
                                  host_directory=hostdir,
                                  remove=remove,
                                  backup=backup,
                                  localname=locname,
                                  local_dir=locdir,
                                  jobid=jobid,
                                  envcmds=envcmd,
                                  verbosity=verbosity)
        except Exception as e:
            print(f"**Encountered Error while putting back to origin machine: {e}")
#            full traceback
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)
            back = False
        if back:
            print("Successfully transferred back the job")
        else:
            print("Couldn't transfer back the job, see the full traceback")

    world.barrier()  # wait for the master to finish.


def savetodb(job, *args, **kwargs):

#    if world.rank == 0: # write only with master
        # with asedbconnect("out.db") as mydb:
        #     mydb.write(job.atoms)
    # parprint("Successfully written into the database")
        job._savetodb(*args, **kwargs)


def main():

    try:
        job = pickle.load(paropen(filename, "rb"))
    except Exception as e:
        print(f"Encountered exception {e} while unpickling the file 'job.p'")
        raise e

    initialize_calc(job)

    try:
        conv = job.run()
        if job.run_type == "static":
            job.converged = job.atoms.calc.scf.converged  # check for scf cyclc
        elif job.run_type == "relax":
            job.converged = conv

    except Exception as e:
        parprint(f"***Encountered Exception {e}")
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_exception(exc_type,exc_value, exc_traceback)
        job.converged = False
        # write the rest of outputs anyway - do not raise the exception

    # save calculator_first
    if savecalc and not job.inputs["calc_args"].get("attach", None):
        try:
            job.atoms.calc.write("calc.gpw")
        except Exception as e:
            parprint(f"***Encountered Exception {e}")
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)

    # write_outputs
    # save db

    if savedb and job.converged: # only save the out.db if the job is converged..
        savetodb(job=job)

    scf = job.atoms.calc.scf
    job.scf = scf
    job.atoms.calc = None
    if job.newrunscheme:
        if hasattr(job.newrunscheme, "calculator"):
            job.newrunscheme.calculator = None

        if hasattr(job.newrunscheme, "atoms"):
            if hasattr(job.newrunscheme.atoms, "calc"):
                if job.newrunscheme.atoms.calc:
                    job.newrunscheme.atoms.calc = None

    with paropen(job.defaults_files["output"], mode="wb") as f: # write only the master
        pickle.dump(job, f)

    if job.autofetch:
        autosendback(remove=remove, backup=backup, envcmd=envcmd)

#         now go back to origin_host
#         TODO: Do this in a better way, It's a work around

if __name__ == "__main__":

    args = commandlineargs()
    filename = args.filename
    savedb = args.savedb in ["True", True]
    savecalc = args.savecalc in ["True", True]
    backup = args.backup in ["True", True]
    remove = args.remove in ["True", True]
    envcmd = None if args.envcmd == "None" else args.envcmd

    try:
        envcmd = envcmd.strip("'") # because we added commas #TODO: better way is to save it in a list or
    # put this as part of the job attributes
    except AttributeError:
        pass

    parprint("Inputs parameters:")
    parprint("filename:{0}\nsave database: {1}, save calculator {2}\nbackup {3}, remove {4}, envcommands {5}".format(filename,
                                                                                                    savedb,
                                                                                                    savecalc,
                                                                                                    backup,
                                                                                                    remove,
                                                                                                    envcmd))

    main()