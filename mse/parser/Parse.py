from collections import OrderedDict
import re
import json
import warnings


def list_dict(input=list):

    """Combines the elements of given input (list) containing dictionaries that spans over two/or more elements
    input: list"""

    outdict = OrderedDict()
    for c,line in enumerate(input):

        line = line.split(":", maxsplit=1)
        key = line[0].strip()
        value = line[-1].strip()

        if value.startswith("{"):
            for secondline in input[c+1::]:
                secondline = secondline.strip()
                value = value.strip() + secondline.strip()
                if secondline.endswith("}"):
                    break
            outdict[key] = value

        elif value.endswith("}") or value.endswith(","):
            continue
        else:
            outdict[key] = value
    return outdict


def parse_dict(input=str):
    out = {}
    input = input.strip()
    if input.startswith("{") and input.endswith("}") and ":" in input:
        input = input.strip("{},")
        input = input.split(",") # This will give us key-value pairs separated
        for k in input:
            k = k.split(":")
            key = k[0]
            value = k[-1]
            if value.startswith("{") and value.endswith("}") and ":" in value:
                value = parse_dict(value)
            else:
                try:
                    value = int(value)
                except ValueError:
                    pass
                try:
                    value = float(value)
                except ValueError:
                    pass
                try:
                    value = list(value)
                except ValueError:
                    pass
                out[key] = value
                return out
    else:
        return input # it is not a dictionary output without any change!


def parse_type_json(input=str, ):
    """Uses Json function to convert a '{key:value}' string into  a dictionary
    :param input: (str) dictionary in the form of string= '{key1: value, key2: value ... }'

    :return: dict
    :"""

    input_org = input

    try:
        out = json.loads(input)

    except AttributeError:
        return input_org
    except json.JSONDecodeError:
        input = re.sub(r'([A-Za-z\-\_\@\*]+\-?[A-Za-z]+)', r'"\1"', input)
        try:

            out = json.loads(input) # will fail for tuple

        except json.JSONDecodeError:
            warnings.warn("Failed to parse the given {} into respective python type".format(input_org))
            return input_org

    return out