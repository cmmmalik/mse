from ase.atoms import Atoms as aseatoms
from ase.db import row


from mse.analysis.aserows import get_calc_parameters, get_args
from mse.workflows.relaxation_scheme import Workflow_relaxation
from mse.workflows.base import Baseworkflow
from mse.utilities import print_nice

import warnings


def set_workflow(c,
                 encut,
                 row=None,
                 atoms=None,
                 dry_run=True,
                 modeargs: dict = None,
                 calcargs: dict = None,
                 check: bool = True,
                 clean_dir: bool = False,
                 symmetrize: bool = False,
                 hpckwargs: dict = None,
                 relaxargs: dict = None):
    warnings.warn("The usage of this function is depracated,  use 'set_workflow_relaxation' instead ", DeprecationWarning)
    warnings.warn("only GPAW is supported")

    if not modeargs:
        modeargs = dict()
    if not calcargs:
        calcargs = dict()
    if not relaxargs:
        relaxargs = dict()

    if row:
        calcs = get_calc_parameters(row)
        calcinps, modeinps = get_args(calcs)
    else:
        calcinps = dict()
        modeinps = dict()

    modeinps.pop("encut", None)
    modeinps["encut"] = encut
    modeinps.update(modeargs)
    calcinps.update(calcargs)

    print(calcinps, modeinps)

    if not atoms and row:
        atoms = row.toatoms(False)
    elif not atoms and not row:
        raise RuntimeError("Either atoms or row must be provided")

    atoms.calc = None

    if not dry_run:
        try:
            wrkf = Workflow_relaxation(atoms=atoms,
                                       working_directory=c,
                                       dircheck=check)

            if symmetrize:
                wrkf.refine_symmetry(symprec=0.001)

            wrkf.initialize_job(name=c,
                                modeinps=modeinps,
                                calcinps=calcinps,
                                relaxinps=relaxargs, )
            if clean_dir:
                wrkf.job.reset()

            print(wrkf.job.relax_inputs)

            wrkf.make_ready(autofetch=True)
            if not hpckwargs:
                hpckwargs = dict()
            server = hpckwargs.pop("server", "lcluster13")
            walltime = hpckwargs.pop("walltime", "12:00:00")
            save_calc = hpckwargs.pop("save_calc", False)
            save_db = hpckwargs.pop("save_db", True)
            remove = hpckwargs.pop("remove", False)
            backup = hpckwargs.pop("backup", False)
            envcmd = "'{}'".format(hpckwargs.pop("envcmd", 'conda activate'))
            wrkf.hpc_setup(server=server,
                           walltime=walltime,
                           save_calc=save_calc,
                           save_db=save_db,
                           remove=remove,
                           backup=backup,
                           envcmd=envcmd, **hpckwargs)
        except Exception as error:
            warnings.warn(f"Got exception:\n{error}")
            try:
                return wrkf
            except NameError:
                raise
        return wrkf


def set_workflow_relaxation(composition: str,
                            working_directory: str,
                            encut: int or float = None,
                            calculator_type: "gpaw" or "vasp" = "gpaw",
                            run_type: "static" or "relax" = "relax",
                            row: row = None,
                            atoms: aseatoms = None,
                            dry_run: bool = True,
                            modeargs: dict = None,
                            calcargs: dict = None,
                            dircheck: bool = True,
                            clean_dir: bool = False,
                            symmetrize: bool = False,
                            symprec: float or int = 0.01,
                            jobargs:dict = None,
                            hpckwargs: dict = None,
                            relaxargs: dict = None,
                            verbosity: int = 1,
                            **extraargs):

        # assert atoms or row

        modeargs = dict() if not modeargs else modeargs
        calcargs = dict() if not calcargs else calcargs
        relaxargs = dict() if not relaxargs else relaxargs
        jobargs = dict() if not jobargs else jobargs

        _argvalues = locals()
        args = ["composition", "working_directory", "encut", "calculator_type", "run_type", "dry_run","dircheck",
                "clean_dir", "symmetrize", "symprec", "verbosity", "jobargs" ]
        _argdct = {k:_argvalues[k] for k in args}


        if row:
            if row.calculator == "gpaw":
                calcs = get_calc_parameters(row)
                calcinps, modeinps = get_args(calcs)
            else:
                raise NotImplementedError("Under developement, "
                                          "don't know how to handle(extract calculator parameters from) the rows")
        else:
            calcinps = dict()
            modeinps = dict()

        modeinps.update(modeargs)
        calcinps.update(calcargs)

        if verbosity >= 1:
            lsep = "\n\t"
            print("calcinps:{0}{1}".format(lsep, print_nice(calcinps, separater=lsep)))
            print("modeinps:{0}{1}".format(lsep, print_nice(modeinps, separater=lsep)))
            print("relaxinps:{0}{1}".format(lsep, print_nice(relaxargs, separater=lsep)))
            print("hpcinps:{0}{1}".format(lsep, print_nice(hpckwargs, separater=lsep)))
            if extraargs:
                print("Extra args:{0}{1}".format(lsep, print_nice(extraargs, separater=lsep)))

            print("Other input paramters:{0}{1}".format(lsep, print_nice(_argdct, separater=lsep)))

        if not atoms and row:
            atoms = row.toatoms(False)

        if dry_run:
            return

        # atoms.calc = None # remove old calculator from atoms
        assert jobargs.get("restart") or atoms # if restart that means atoms will be read fomr calc.gpw

        if calculator_type == "gpaw":
            if encut:
                modeinps["encut"] = encut

            wrkf = Workflow_relaxation(atoms=atoms,
                                        working_directory=working_directory,
                                       reltype=relaxargs.get("reltype", ["ions","cell", "full"]),
                                        dircheck=dircheck)

        elif calculator_type == "vasp":  # for vasp, we do not need workflow_relaxation,
                # as vasp handles relaxation on its own
            if encut:
                calcinps["encut"] = encut
            wrkf = Baseworkflow(atoms=atoms,
                                calculator_type=calculator_type,
                                working_directory=working_directory,
                                dircheck=dircheck)
            if verbosity >= 1:
                sst = []
                disp = []
                if modeinps:
                    sst.append(["'modeinps'"])
                    disp.append(["modeinps:{}".format(modeinps)])
                if relaxargs:
                    sst.append(["'relaxargs'"])
                    disp.append(["'relaxargs':{}".format(relaxargs)])
                sst = " & ".join(sst)
                disp = "\n".join(disp)

                warnings.warn(sst + "will not be used for vasp\n" + disp, RuntimeWarning)

                print("For vasp, they belong to 'calcinps'")
                print("Cleaning...")
                print("Calculation_parameters:\n{}".format(calcinps))
            #cleaning modeinps
            modeinps = None
            relaxargs = None

        else:
            raise ValueError("Currenly only two types of calculator are supported: 'vasp' and 'gpaw',"
                             "instead {} was chosen".format(calculator_type))

        if symmetrize:
            wrkf.refine_symmetry(symprec=symprec)

        wrkf.initialize_job(name=composition,
                            run_type=run_type,
                            modeinps=modeinps,
                            calcinps=calcinps,
                            relaxinps=relaxargs,
                            **jobargs)
        if clean_dir:
            wrkf.job.reset()

        wrkf.make_ready(autofetch=extraargs.get("autofetch", True))

        if not hpckwargs:
            hpckwargs = dict()

        server = hpckwargs.pop("server", "lcluster13")
        walltime = hpckwargs.pop("walltime", "12:00:00")
        save_calc = hpckwargs.pop("save_calc", False)
        save_db = hpckwargs.pop("save_db", True)
        remove = hpckwargs.pop("remove", False)
        backup = hpckwargs.pop("backup", False)

        wrkf.hpc_setup(server=server,
                       walltime=walltime,
                       save_calc=save_calc,
                       save_db=save_db,
                       remove=remove,
                       backup=backup,
                       **hpckwargs)
        return wrkf
