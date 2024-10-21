"""
Utils for snakemake wrapper light-weight system
===============================================

"""

from typing import Union, Callable, Optional
from snakemake.io import Params
from snakemake.utils import update_config


__author__ = "Fati Chen & Julie Ripoll"
__created__ = "2020-05"
__license__ = "CeCILL"
__version__ = "0.1.0"
__vscript__ = "0.1.0"


class _StubParam(str):
    """ 
    Represents a ghost parameter, will return False and is equal to None is represented as a snakemake object eg: {params.example}
    
    """
    __stubbed = True

    def __eq__(self, other):
        return None == other

    def __nonzero__(self):
        return False

    def __bool__(self):
        return False


class _StubParams(Params):
    """ 
    Represents an empty Params object, will return False or will delegate resolution to snakemake
    
    """
    __stubbed = True

    def __init__(self, prefix):
        self.__prefix = prefix

    def __stub_param(self, key):
        return _StubParam("{"+self.__prefix+"."+key+"}")

    def __getitem__(self, key):
        return self.__stub_param(key)

    def __getattr__(self, key):
        return self.__stub_param(key)

    def __str__(self):
        return "{"+self.__prefix+"}"

    def get(self, key, default=""):
        return default


STUB_SNAKEMAKE_VARS = {
    "input": _StubParams("input"),
    "log": _StubParams("log"),
    "output": _StubParams("output"),
    "params": _StubParams("params"),
    "resources": _StubParams("resources"),
    "rule": _StubParams("rule"),
    "threads": _StubParams("threads"),
    "wildcards": _StubParams("wildcards")
}


def inject(fn, snakemake_vars):
    """ Injects snakemake variables into the function.

    If snakemake vars cannot be resolved, injects a ghost object which impersonate snakemake variables, always returning false
    
    """
    if 'rule' in snakemake_vars:
        shell_command = fn(**snakemake_vars)
        print(shell_command)
        return shell_command
    else:
        return "WARNING : snakemake values are not currently available, can result in inconsistencies in the displayed shell command.\n"+fn(**STUB_SNAKEMAKE_VARS)


def merge_dict(default, overwrite):
    """Merges default and overwrite into an new object.

    .. deprecated::
        `merge_dict` has too many limitations and will be removed in the next iterations
    
    """
    if isinstance(overwrite, _StubParams):
        return Params(fromdict=default)

    config = {}
    update_config(config, default)
    update_config(config, overwrite)
    return Params(fromdict=config)


def merge_strings(*args, join_string=" "):  # old join " \\\n"
    """
    .. deprecated::
        `merge_strings` is now deprecated and will be removed in the next iterations, use
        `join_str` instead
    
    """
    return join_str(*args, joiner=join_string)


def join_str(*args, joiner=" "):  # old join " \\\n" and test results 'a \\\nb \\\nc'
    r""" Joins argument strings using joiner attribute.

    Will ignore Falsy strings

    Parameters
    ----------
        *args: str[]
            strings to join together
        joiner: str
            string used to join `*args`, default allows newline for cmds (default: '\\\\n ')

    Returns
    -------
    the joined string

    Examples
    --------

        >>> join_str('a', 'b', 'c')
        'a b c'

        >>> join_str('a', 'b', 'c', joiner=', ')
        'a, b, c'
    
    """
    return joiner.join(filter(None, args))


def val_mappy(value: str, matcher: Union[str, dict, tuple, list, Callable[..., Union[str, None, bool]], None] = None, default: Optional[str] = None):
    """Will try to match value with the matcher.

    If `default` is not provided and val_mappy cannot match the variable. it will raise an error.

    Note
    ----
    During DAG building phase, using val_mappy results 

    Parameters
    ----------
        value: str
            the value to resolve using the matcher

        matcher: str, dict, tuple, list, callable 
            holds the values to compare to.
            Depending on the type of the matcher, it will behave differently:

            str: will check if `value == matcher`, if so returns `value`
                returns `default` or raises ValueError otherwise

            dict: if `value` is present as a key, will return the corresponding `value`
                returns `default` or raises ValueError otherwise

            tuple, list: checks if `value` is present in the tuple/list
                returns `default` or raises ValueError otherwise

            callable:  (`value`: str, `default`: str) -> bool, str, None
                if the result is not None, val_mappy returns the result
                returns `default` or raises ValueError otherwise

        default: str, None
            defined, it will be returned if value cannot resolve using matcher.
            if None, mappy will raise ValueError if value is None

    Examples
    --------

        >>> params = Params(fromdict={"param1": "value 1", "param2": "value2"})

    You can use it to require the existance of a parameter

        >>> val_mappy(params.param1)
        'value 1'

        >>> val_mappy(params.get('param5')) is None
        True

    Exact match

        >>> val_mappy(params.param1, "value 1")
        'value 1'

        >>> val_mappy(params.param1, "does not exist")
        Traceback (most recent call last):
        ...
        ValueError: ("'value 1' does not match given value and no default value have been defined", 'value 1', 'does not exist')

    Using a tuple of allowed values

        >>> val_mappy(params.param1, ("value 1", "other allowed value"))
        'value 1'

        >>> val_mappy(params.param1, ("value11", "other allowed value"))
        Traceback (most recent call last):
        ...
        ValueError: ("'value 1' does not match given values and no default value have been defined", 'value 1', ('value11', 'other allowed value'))

    Using a list of allowed values

        >>> val_mappy(params.param1, ["value 1", "other value"])
        'value 1'

        >>> val_mappy(params.param1, ["value11", "other value"])
        Traceback (most recent call last):
        ...
        ValueError: ("'value 1' does not match given values and no default value have been defined", 'value 1', ['value11', 'other value'])

    Using a dictionary for existing values and mapping

        >>> val_mappy(params.param1, {'value 1': 'another object', 'other accepted': 'Hello World'})
        'another object'

        >>> val_mappy(params.param1, {'value11': 'will raise an error'})
        Traceback (most recent call last):
        ...
        ValueError: ("'value 1' does not match given values and no default value have been defined", 'value 1', {'value11': 'will raise an error'})

    Using a function

        >>> val_mappy(params.param1, lambda needle, default : needle+" augmented")
        'value 1 augmented'

        >>> def validation(needle, default):
        ...     return needle if needle.startswith('value') else None
        >>> val_mappy(params.param1, validation)
        'value 1'

    See Also
    --------

    mappy: solving snakemake arguments
    
    """
    try:
        value.__StubParam__stubbed
    except AttributeError:
        pass
    else:
        mapper_str = f"{{{matcher!r}}}" if isinstance(
            matcher, dict) else f'{matcher!r}'
        default_str = "" if default == None else f', default={default!r}'
        return f"<val_mappy({value}, {mapper_str}, {default_str})>"

    if isinstance(matcher, str):
        if matcher == value:
            return value
        else:
            if default == None:
                error_msg = f"{value!r} does not match given value and no default value have been defined"
                raise ValueError(error_msg, value, matcher)
            return default

    elif isinstance(matcher, dict):
        try:
            return matcher[value]
        except KeyError:
            pass

        if default == None:
            error_msg = f"{value!r} does not match given values and no default value have been defined"
            raise ValueError(error_msg, value, matcher)
        return default

    elif isinstance(matcher, (tuple, list)):
        if value in matcher:
            return value
        else:
            if default == None:
                error_msg = f"{value!r} does not match given values and no default value have been defined"
                raise ValueError(error_msg, value, matcher)

            return default

    elif callable(matcher):
        result = matcher(value, default)
        if result != None:
            return result

        if default != None:
            return default

        error_msg = f"passing {value!r} through {matcher!r} returned None and no default value have been defined"
        raise ValueError(error_msg, value, matcher)

    return value


def mappy(holder, variable: str, matcher: Union[str, dict, tuple, list, Callable[..., Union[str, None, bool]], None] = None, default=None):
    """ Will try to retrieve the `variable` attribute of `holder` using `holder.get(variable, default)` and match it with the matcher.

    If default` value is not provided and mappy cannot find/match the variable. it will raise an error.

    Parameters
    ----------
        holder : snakemake input, output, params ...
            the object containing the variable

        variable: str
            name of the variable, the value (`needle`) of `holder.get(variable)` will be retrieved .
            if holder does not contain `variable`, will use `default` as the needle.

        matcher: str, dict, tuple, list, callable 
            holds the values to compare to.
            Depending on the type of the matcher, it will behave differently:

            str: will check if `needle == matcher`, if so returns needle
                returns `default` or raises ValueError otherwise

            dict: if needle is present as a key, will return the corresponding value
                returns `default` or raises ValueError otherwise

            tuple, list: checks if `needle` is present in the tuple/list
                returns `default` or raises ValueError otherwise

            callable:  (needle: str, default: str) -> bool, str, None
                in the case that holder does not contain variable, needle and default will hold the same value.
                if the result of the function is not None, mappy will return the result
                returns `default` or raises ValueError otherwise
        default: str, None
            defined, it will be returned in case that holder does not contain variable, or if it cannot resolve using matcher.
            if None, mappy will raise ValueError if needle is None

    Returns
    -------
    Any
        the resolved value using matcher, 
        if default is defined and mappy could not resolve the value, will return default
        raises ValueError if unable to resolve and default is None

    Examples
    --------

        >>> params = Params(fromdict={"param1": "value 1", "param2": "value2"})

    You can use it to require the existance of a parameter

        >>> mappy(params, "param1")
        'value 1'

        >>> mappy(params, "param5")
        Traceback (most recent call last):
        ...
        ValueError: the value for 'param5' is None, you can avoid this error by defining a default value

    Exact match

        >>> mappy(params, "param1", "value 1")
        'value 1'

        >>> mappy(params, "param1", "does not exist")
        Traceback (most recent call last):
        ...
        ValueError: ("'value 1' does not match given value and no default value have been defined", 'value 1', 'does not exist')

    Using a tuple of allowed values

        >>> mappy(params, "param1", ("value 1", "other allowed value"))
        'value 1'

        >>> mappy(params, "param1", ("value11", "other allowed value"))
        Traceback (most recent call last):
        ...
        ValueError: ("'value 1' does not match given values and no default value have been defined", 'value 1', ('value11', 'other allowed value'))

    Using a list of allowed values

        >>> mappy(params, "param1", ["value 1", "other value"])
        'value 1'

        >>> mappy(params, "param1", ["value11", "other value"])
        Traceback (most recent call last):
        ...
        ValueError: ("'value 1' does not match given values and no default value have been defined", 'value 1', ['value11', 'other value'])

    Using a dictionary for existing values and mapping

        >>> mappy(params, 'param1', {'value 1': 'another object', 'other accepted': 'Hello World'})
        'another object'

        >>> mappy(params, 'param1', {'value11': 'will raise an error'})
        Traceback (most recent call last):
        ...
        ValueError: ("'value 1' does not match given values and no default value have been defined", 'value 1', {'value11': 'will raise an error'})

    Using a function

        >>> mappy(params, 'param1', lambda needle, default : needle+" augmented")
        'value 1 augmented'

        >>> def validation(needle, default):
        ...     return needle if needle.startswith('value') else None
        >>> mappy(params, 'param1', validation)
        'value 1'

    See Also
    --------

    val_mappy: directly solving value against the matcher
    
    """
    try:
        holder.__StubParams__stubbed
    except AttributeError:
        pass
    else:
        matcher_protect = f"{matcher!r}".replace("{", "{{").replace("}", "}}")
        mapper_str = f"{matcher_protect}" if isinstance(
            matcher, dict) else f'{matcher!r}'
        default_str = "" if default == None else f', default={default!r}'
        return f"<mappy({holder[variable]!r}, {mapper_str}{default_str})>"

    value = holder.get(variable, default)

    if value == None:  # means that value and default are none
        error_msg = f"the value for {variable!r} is None, you can avoid this error by defining a default value"
        raise ValueError(error_msg)

    return val_mappy(value, matcher, default=default)


def optional(holder, variable, string, default_value=None):
    """
    Checks if holder has variable, will return a formatted `string` with the corresponding value injected

    Parameters
    ----------
        holder : snakemake input, output, params ...
            the object containing the variable
        string : str
            the value to format and return
        variable : str
            the name of the variable
        default_value : str, optional
            default value when variable is not set (default : None)

    Returns
    -------
    the formatted string, 
    if it is called with stubbed parameters, will format the string with `<variable>`

    Examples
    --------
        >>> input = Params(fromdict={"annot": "example.gff"})
        >>> optional(input, "annot", "--annotation {}")
        '--annotation example.gff'

        >>> optional(input, "unexistant", "--unexistant {}")
        ''

        During snakemake printshellcmds:

        >>> stub_input = _StubParams("input")

        >>> optional(stub_input, "annot", "--annotation {}")
        '[--annotation <{{input.annot}}, default=None>]'
    
    """
    try:
        holder._StubParams__stubbed
    except AttributeError:
        pass
    else:
        return "[{}]".format(string.format(f"<{{{holder[variable]}}}, default={default_value!r}>"))

    value = holder.get(variable, default_value)
    return string.format(value) if value else ""


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
