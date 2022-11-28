import argparse
import os
import sys
import traceback
import time
from mse.servertools.ssh import ShellHandler
import warnings

from HPCtools.IO import Logger


def wrap_path(wraper="gpaw_wrap", host=None):

    """
    @param wraper: str, the name of wrapper script
    @param host: str, hostname of the server, can be hostnames from config file
    @return:  str, path where the gpaw_wrap script is located on the HPC server.
    """
    warnings.warn()
    terminal = ShellHandler(host=host)
    #load modules
    _, out, err = terminal.clean_execute("module load intel/2018\nmodule load intelmpi\n"
                                         "module load python\n"
                                         "module load gpaw") # need to load these modules # otherwise we will
    #run into an error.

    if err:
        raise RuntimeError(f"Encountered error: {err}")
    _, out, err = terminal.clean_execute(f"python -c 'from mse.wrapper import {wraper}\n" 
                                             "print(gpaw_wrap.__file__)'")
    if err:
        raise RuntimeError(f"Encountered an error: {err}")

    path = out[-1]
    path = path.strip()
    return path

def get_timeout(directory="."):
    size = os.path.getsize(directory)
    size = size/1000 #kByte
    if size <= 1e2:
        return size + 20
    return 15*60

def putback_origin(jobname:str,
                   host: str,
                   host_directory: str,
                   local_dir: str,
                   localname: str,
                   remove: bool,
                   backup: bool,
                   envcmds: str = None,
                   jobid: int = None,
                   timeout = 250,
                   verbosity: int = 0):

#   TODO: Better way would be to have the origin of the job script, retrieving back upon the job completion/finish,
#   this is a work around at the moment.

    """

    @param host: must be present in the ssh config file and passwordless login must be enabled to the ending machine.
    @param host_directory:
    @param jobname:
    @param remote:
    @envcmds: (str) extra environmental(bash) commands that must be invoked before exceuting the actual script. e.g.
    'conda activate' etc.
    @return:
    """

    terminal = ShellHandler(host=host) # Uses Paramiko terminal instead of naive ssh.

    _, out, err = terminal.clean_execute("cd %s" % host_directory)
    if envcmds:
        if verbosity >= 2:
            print("Debug: envcmd being executed\n: {}".format(envcmds))
        _, out, err = terminal.clean_execute(envcmds)
# TODO: Put this in an argument in extras, as environmental setups

    #cmd += [f"python -c 'from mse.Jobs.job import HPCjob\n" + \
     #     f"hpc = HPCjob(jobname=None, working_directory=\"{host_directory}\")\n" + \
     #     f"print(\"Fetching....\")\n" + \
     #     f"hpc.fetch_hpc()'"]

    cmd = "python -u -m mse.servertools.Slurm"
    cmd += " -p {}".format(host_directory)
    cmd += " -b {}".format(backup)
    cmd += " -r {}".format(remove)

    if jobname is not None:
        cmd += " -j {}".format(jobname)

    if jobid is not None:
        cmd += " -ji {}".format(jobid)

    cmd += " -H {}".format(localname)
    cmd += " -Hd {}".format(local_dir)
    cmd += " -v {}".format(verbosity)
    cmd += " &"

    if verbosity >= 1:
        print("Executing command on the server for fetching back\n{}".format(cmd))

    terminal.execute(cmd)
    exit_code = 1 # Assume an error until unless proved otherwise

    # set timeout
    if not timeout:
        timeout = get_timeout(directory=".")

    terminal.channel.settimeout(timeout)

    for line in terminal.stdout: # There is a chance that we can get stuck here! # Add timeout limit
        print(line, flush=True)
        if "Exit code" in line:
            exit_code = int(line.rsplit(maxsplit=1)[-1])
            break
        else:
            time.sleep(1) # wait for the output

    terminal.close()

    if exit_code == 1:
        print(f"***Encountered Exception")
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_exception(exc_type,exc_value, exc_traceback)
        raise RuntimeError("Encountered an Error:{}".format(exit_code))

    return True


def commandlineargs():

    parser = argparse.ArgumentParser()
    parser.add_argument("--path", "-p",
                        help="Path of working_directory of the job",
                        default=os.getcwd(),
                        type=str,
                        dest="working_directory")

    parser.add_argument("--remove", "-r",
                        help="set True(or False) to remove(or not ) the folder",
                        default=False,
                        type=str,
                        dest="remove")

    parser.add_argument("--backup", "-b",
                        help="Set to True(or False) to backup(or not) the local folder contents",
                        default=False,
                        type=str,
                        dest="backup")

    parser.add_argument("--jobname", "-j",
                        help="Enter the name of the job(not used or required"
                             "in practice, only for all inclusivity purposes)",
                        dest="jobname")

    parser.add_argument("--jobid", "-ji",
                        help="Job id(HPC) of the job to update the database in the dummy host machine (=from where the job )",
                        dest="jobid",
                        type=int)

    parser.add_argument("--hostname", "-H",
                        help="Actual Host name (on which the job was running)",
                        dest="hostname")

    parser.add_argument("--hpcdir", "-Hd",
                        help="Actual host directory",
                        dest="hpcdir")

    parser.add_argument("--verbosity", "-v",
                        help="Sets Verbosity level (0=default),1,2 -- useful for debugging the code",
                        default=0,
                        type=int,
                        dest="verbosity")

    args = parser.parse_args()
    return args


def main():

    from mse.Jobs.job import HPCjob
    sys.stdout.write("Fetching back")
    sys.stdout.flush()

    if not jobid or not hpcdir or not hostname:

        warnings.warn("Certain parameters were not provided,"
                      " will fall back to reading the information from the working_directory", RuntimeWarning)
        try:
            warnings.filterwarnings("error")
            hpc = HPCjob.read_hpc(jobname=jobname,
                                working_directory=working_directory,
                                jobid=jobid,
                                status="Running",
                                localdir=working_directory)

        except UserWarning as e:
            warnings.filterwarnings("default") # TODO: Need to test this fragment of code
            raise UserWarning("iid/jobid was not provided to look into the database, job cannot be fetched back")

    else:
        hpc = HPCjob.from_server_hpc(jobname=jobname,
                                    working_directory=working_directory,
                                    hostname=hostname,
                                    jobid=jobid,
                                    hpcdir=hpcdir,
                                    localdir=working_directory,
                                    )

    print(hpc)

    try:
        print("Options:\nremove set to {0}, backup set to {1}".format(remove, backup))
        hpc.fetch_hpc(status="Completed",
                      remove=remove,
                      backup=backup)
        # if we get here , we assumed that job part is done

    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_traceback)
        sys.stdout.write(f"**Encountered Excepion {e}")
        sys.stdout.write("Exit code 1\n")
        raise e

    sys.stdout.write("Exit code 0\n")


if __name__ == "__main__":

    args = commandlineargs()
    working_directory = args.working_directory
    remove = args.remove in ["True", True]
    backup = args.backup in ["True", True]
    jobname = None if args.jobname == "None" else args.jobname
    jobid = None if args.jobid == "None" else args.jobid
    hostname = None if args.hostname == "None" else args.hostname
    hpcdir = None if args.hpcdir == "None" else args.hpcdir
    verbosity = args.verbosity

    main()


