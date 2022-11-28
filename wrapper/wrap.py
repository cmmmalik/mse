from ase.parallel import world, parprint as aseparprint, paropen
import sys
from functools import partial
import traceback

from mse.servertools.Slurm import putback_origin
from HPCtools.hpc_tools3 import HPCMain

parprint = partial(aseparprint, flush=True) # flush the buffers
print = partial(print, flush=True)


def savetodb(job):
    if world.rank == 0: # write only with master
        job._savetodb()


def autosendback(remove: bool,
                 backup: bool,
                 envcmd: str,
                 verbosity: int = 0):

# Host refers to non-local machine, present on the other end of ssh connection.
# Do not raise any exception, let the job pass through to completion.

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