import re


def remove_1(st):
    """Removes the character '1' in the input string"""
    matched = re.findall("[A-Z][a-z]?1(?=\D)", st) + re.findall("[A-Z][a-z]?1$", st)

    for i in matched:
        loc = st.find(i)
        st = st.replace(i, i[:-1])
    return st


def float_int(number):
    """returns int-type string representation of the floating number that contains only tailling zeros after the decimals
     Useful for fine string printing"""

    matched = re.search("(?<=\d\.)0+(?!\d+)", str(number))
    if matched:
        number = int(number)
    return number