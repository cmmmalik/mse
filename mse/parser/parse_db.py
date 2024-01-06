
def keys_parsedb(data:dict):
    keys = {k: (str(value) if isinstance(value, dict) or
                              isinstance(value, list) or
                              isinstance(value, tuple) else value) for k, value in data.items()}
    return keys