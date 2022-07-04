#!/usr/bin/env python3

from argparse import ArgumentParser

from mse.Jobs.job import HPCjob

def args():
    parser = ArgumentParser()
    parser.add_argument("-r", "--read", help="supply if to read the hpc from the directory", dest="read", nargs="?", const=True, default=True)
    parser.add_argument("-w", "--workdir", help="Working directory", action="store", dest="wrk", type=str, default=None)
    parser.add_argument("-v", "--verbosity", help="verbosity", action="count", dest="verbosity", default=0)
    parser.add_argument("-f", "--fetch", help="Fetch the job", action="store", dest="fetch", default=False, type=bool)
    parser.add_argument("-e", "--extra-args", help="Extra arguments that will be provided to a given operation", action="store",
                        nargs="*", default=[], dest="kwargs")

    args = parser.parse_args()
    return args


def parse_extraarguments(arg):

    kwdict = {}
    for keyvalue in arg:
        k, v = keyvalue.split("=")
        if "=" in k or "=" in v:
            raise ValueError("INvalid command line argument provided")
        elif "True" in v or "true" in v:
            v = True
        elif v.isalpha():
             pass
        else:
            try:
                v = int(v)
            except ValueError:
                v = float(v)
        kwdict[k] = v
        return kwdict


def main():
    ags = args()
    working_directory = ags.wrk
    read = ags.read
    verbosity = ags.verbosity
    fetch = ags.fetch
    kwargs = parse_extraarguments(ags.kwargs)
    if kwargs is None:
        kwargs = {}

    if read is True:
        print("Reading the hpc ....")
        hpc = HPCjob.read_hpc(working_directory=working_directory)
        if verbosity >= 1:
            print(hpc)

        if fetch is True:
            print("Fetching ...")
            hpc.fetch_hpc(**kwargs)


if __name__ == "__main__":
    main()