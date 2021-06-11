from ase.parallel import paropen
import argparse
import sys
import pickle
import traceback

from Simtools.wrapper.wrap import savetodb, parprint, print, autosendback


def generate_runcmd(inputfile,
                    save_calc: bool = True,
                    save_db: bool = False,
                    backup: bool = False,
                    remove: bool = True,
                    envcmd: bool = None):

    loc = "Simtools.wrapper.vasp_wrap"
    runcmd = "python -m {}".format(loc)
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


def main():
    try:
        job = pickle.load(paropen(filename, "rb"))
    except Exception as e:
        print(f"Encountered exception {e} while unpickling the file 'job.p'")
        raise e

    job.initialize_calc()

    try:
        out = job.run()
        conv = out is not None
        job._converged = conv
    except Exception as e:
        parprint(f"***Encountered Exception {e}")
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_exception(exc_type,exc_value, exc_traceback)
        job._converged = False

    if savedb:
        savetodb(job)

    if savecalc:
        try:
            job.atoms.calc.write_json("calc.json")
        except Exception as e:
            parprint(f"***Encountered Exception {e}")
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)

    job.atoms.calc = None
    with paropen(job.default_files["output"], mode="wb") as f:
        pickle.dump(job, f)

    if job.autofetch:
        autosendback(remove=remove,
                     backup=backup,
                     envcmd=envcmd,
                     )

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