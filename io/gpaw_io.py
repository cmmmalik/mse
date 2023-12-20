import warnings
import re

from mse.parser.Parse import list_dict,  parse_type_json

from copy import deepcopy
import json


class FormatException(Exception):
    pass


class Readparameters:

    """Read initial calculation parameters from output.txt file"""

    def __init__(self, file, format="gpaw-out", index=None):
        if format != "gpaw-out":
            raise FormatException("Currently, gpaw-out text file is supported")

        self.file = file
        self.inputs = None
        self.name = "gpaw"
        self.version = None
        self.subpackages = None
        self.parameters = {}
        self._indexparameters = []
        self._indexinputs = []
        self.kp = None

        with open(file, "r") as txtfile:

            self._grab_version(txtfile)
            self._subpackages_info(txtfile)

            for line in txtfile:
                if line.startswith("Input parameters"):
                    self.inputs = []
                    self._grab_input_parameters(txtfile) # No need to convert it to dictionary at the moment
                    self._indexinputs.append(self.inputs)
                # Output(written) parameters

                elif line.startswith("Occupation numbers:"): #Sometimes, we don't have this line, in case of default,
                    #written after wavefunction info
                    if line.endswith("numbers:"):
                        nxtline = txtfile.__next__()
                    else:
                        nxtline = line

                    nxtline = nxtline.strip().split(":", maxsplit=1)
                    # if self.parameters is None:
                    #     self.parameters = {}
                    # else:
                    #     self._indexparameters.append(self.parameters.copy())

                    self.parameters["occupations"] = nxtline # No need to parse it!

                elif line.startswith(" "*2 + "Fermi-Dirac:"):
                    self.parameters["occupations"] = line.strip()

                #Kpoints
                elif "k-points:" in line and "Parallelization over" not in line:
                    line = line.split(" ", maxsplit=1)
                    try:
                        k_no = int(line[0].strip())
                    except ValueError:
                        k_no = None
                    line = line[-1].split(":", maxsplit=1)
                    self.parameters[line[0]] = line[-1].strip()
                elif "Wave functions:" in line:
                    self.parameters["Wave functions"] = line.split(":", maxsplit=1)[-1].strip()
                    self._wavefunction_info(txtfile)
                elif line.startswith("Eigensolver"):
                    nxtline = txtfile.__next__()
                    self.parameters["eigensolver"] = nxtline.strip()
                elif line.startswith("Hamiltonian:"):
                    self._getxc_info(txtfile)

                elif line.startswith("Fermi level"):
                    self._indexparameters.append(self.parameters.copy())

            if index:
                self.parameters = self._indexparameters[index]
                try:

                    self.inputs = self._indexinputs[index]
                except IndexError:
                    self.inputs = None

            if self.inputs:
                self.inputs = list_dict(self.inputs)
                self._parse_input_parameters()

            self._parse_occupations()
            self._parse_k_points()

    def __repr__(self):
        st = "{self.name:<8}, parameters={self.parameters}," \
             "inputs(txtfile)={self.inputs}," \
             "kp(.py)={self.kp}"

        return st.format(self=self)

    def __str__(self):
        return self.__repr__()

    def _getxc_info(self, txtfile):
        for nxtline in txtfile:
            if nxtline.strip().endswith("Correlation functional"):
                # s = re.search("(?<=\s)[A-Z]+(?=\sExchange)", nxtline)
                s = re.search("(?<=\s)\S*(?=\sExchange)", nxtline)
                if not s:
                    warnings.warn("Couldn't find any match for 'Exchange-Correlation', (debug it)")
                    break
                self.parameters["xc"] = s.group()
                break

    def _wavefunction_info(self, txtfile):

        for counter in range(0, 3):
            nxtline = txtfile.__next__()
            if nxtline.strip().startswith("Cutoff"):
                nxtline = nxtline.split(":", maxsplit=1)
                nxtline = nxtline[-1].split(maxsplit=1)
                self.parameters["encut"] = float(nxtline[0])
            elif nxtline.strip().startswith("Pulay-stresss"):
                nxtline = nxtline.split(":",maxsplit=1)
                nxtline = nxtline[-1].strip().split(" ")
                self.parameters["dedecut"] = float(nxtline[0])
                break

    def _subpackages_info(self, txtfile):
        if self.subpackages is None:
            self.subpackages = {}
        for line in txtfile:
            if line.startswith("User:"):
                line = line.rsplit(":", maxsplit=1)
                self.subpackages[line[0].lower()] = line[-1].strip() #save the first line
                for nextline in txtfile:
                    if nextline == "\n":
                        break
                    if nextline.startswith("_gpaw") or nextline.startswith("  "):
                        continue

                    nextline = nextline.rsplit(":", maxsplit=1)
                    self.subpackages[nextline[0].lower()] = nextline[-1].strip()
    #                if nextline[0].startswith("OMP_NUM"):
    #                    break

                break

    def _grab_version(self, txtfile):

            if not self.version is None:
                warnings.warn("GPAW version is already present", RuntimeWarning)
            for line in txtfile:
                line = line.strip()
                if line.startswith("|__"):
                    line = line.rsplit(maxsplit=1)
                    version = line[-1]
                    if "." in version:
                        self.version = version
                        break

    def _grab_input_parameters(self, openedfile):

        for line in openedfile:
            if line == "\n": # break at empty line
                break
            self.inputs.append(line.strip())

        #self.inputs = ",".join(self.inputs)

    def _parse_input_parameters(self):
        if not isinstance(self.inputs, dict):
            raise RuntimeError("Attribute {} is not an instance of {}".format('self.inputs', dict))
        for k in self.inputs:
            self.inputs[k] = parse_type_json(self.inputs[k])

    def _parse_occupations(self):
        try:
            occu = self.parameters["occupations"]
        except (KeyError, TypeError) as error:
            warnings.warn("attribute {} does not have 'occupations' ".format('self.parameters'))
            return
        out = {}
        if isinstance(occu, str):
            if ":" in occu:
                occu = occu.split(":")
            else:
                warnings.warn(f"Unexpected string '{occu}', unable to parse it, a possible bug!!!")
                return

        out["name"] = occu[0]

        for value in occu[1].split(","):
            if "=" in value:

                value = value.split("=")
                value[1] = re.sub("[a-zA-Z]+?","", value[1])

                try:
                    value[1] = int(value[1])
                except ValueError:
                    value[1] = float(value[1])

                out[value[0].strip()] = value[1]
            else:
                out[value] = value

        self.parameters["occupations"] = out.copy()

    def _parse_k_points(self):
        try:
            kpars = self.parameters["k-points"]
        except (KeyError, TypeError) as error:
            warnings.warn("attribute {} does not have 'k-points'".format('self.parameters'))
            return
        out = {}
        kpts = re.findall("\d+ x \d+ x \d+", kpars)[0]
        kpars = kpars.replace(kpts, "")
        kpts = list(kpts.split("x"))
        out["k-points"] = kpts
        #looking for shift
        if "+" in kpars:
            shift = kpars.split("+")
            kpars = kpars.replace(shift[1], "").replace("+","")
            out["shift"] = shift[1]

        out["name"] = kpars

        self.parameters["k-points"] = out.copy()


def getkdenpy(filename, verbosity: int = 0):
    output = {}
    with open(filename,"r") as f:
        for line in f:
            if line.startswith("kden"):
                line = line.strip()
                line = line.split("=", maxsplit=1)
                try:
                    value = int(line[-1])
                except ValueError:
                    value = float(line[-1])
                output["kden"] = value

            elif line.startswith("kpden") or line.startswith("kpts"):
                line = line.strip()
                line = line.split("=", maxsplit=1)
                value = line[-1]
                if verbosity >=2:
                    print("Debug:\nbefore replacing quotes:\n{}".format(value))
                value = value.replace("\"", "").replace("\'","") # remove old \" around the key values
                if verbosity >= 2:
                    print("Debug:\nbefore parsing \n{}".format(value))
                output["kpden"] = parse_type_json(value)
                # this will come after the first if so the loop can be broken.
                break

    return output


def parse_string_dict(dctstr: str):
    dctstr = re.sub(r"(\')", r'"', dctstr)
    return json.loads(dctstr)

def kpden_gpaw_rows(row):
    kpden = row.calculator_parameters.get("kpts", None)
    if not kpden: # old_row data
        kpden = parse_string_dict(deepcopy(row.kpden))
        kden = row.kden
        kpden["density"] = kden
    else:
        kpden = kpden.copy()
    assert "density" in kpden
    gamma = kpden.pop("gamma", None)
    gamma = True if gamma in ["True", True] else None
    print(kpden, gamma)
    return kpden, gamma