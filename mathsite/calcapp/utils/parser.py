from sympy import *
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication_application,
    convert_xor,
    function_exponentiation,
)
from sympy.physics.units.quantities import Quantity

# English alphabet, both cases
english_letters = [chr(i) for i in range(ord("a"), ord("z") + 1)] + [
    chr(i) for i in range(ord("A"), ord("Z") + 1)
]

# Can add more later
greek_letters = [
    "alpha",
    "beta",
    "lambda",
    "mu",
    "phi",
]

allowed_functions = {
    "sin": sin,
    "cos": cos,
    "tan": tan,
    "asin": asin,
    "acos": acos,
    "atan": atan,
    "sinh": sinh,
    "cosh": cosh,
    "tanh": tanh,
    "exp": exp,
    "ln": log,
    "log": log,
    "sqrt": sqrt,
}

allowed_constants = {
    "pi": pi,
    "e": E,
}

transformations = (
    standard_transformations
    + (implicit_multiplication_application,)
    + (convert_xor,)
    + (function_exponentiation,)
)


def build_local_dict(parameters):
    """
    Build a SymPy local_dict where:
    - parameters (variables): real=True
    - all other alphabetical names: real=True, positive=True
    """
    variables = set(str(p) for p in parameters)

    local = {}

    for name in english_letters + greek_letters:
        if name in variables:
            local[name] = Symbol(name, real=True)
        else:
            local[name] = Symbol(name, real=True, positive=True)

    local.update(allowed_functions)
    local.update(allowed_constants)

    return local


def parse_input(expr_string: str, parameters: list):
    """
    Parse plaintext input into a SymPy expression:
    - variables: real=True
    - constants: real=True, positive=True
    """
    expr_string = expr_string.strip()

    try:
        local_dict = build_local_dict(parameters)

        return parse_expr(
            expr_string,
            transformations=transformations,
            local_dict=local_dict,
            evaluate=True,
        )
    except Exception as e:
        raise ValueError(f"Invalid expression: {expr_string}\nError: {e}")


def parse_bounds(bounds: list, parameters: list):
    """
    Parse a list of two string bounds:
        ["0", "sqrt(2)"]
        ["a", "2*a"]
        ["-pi/2", "pi/2"]

    Returns a list of two SymPy expressions, with correct
    variable/constant assumptions.
    """

    if not isinstance(bounds, list) or len(bounds) != 2:
        raise ValueError("Bounds must be a list of exactly two strings.")

    lower_str, upper_str = bounds

    if not isinstance(lower_str, str) or not isinstance(upper_str, str):
        raise ValueError("Both bounds must be strings.")

    local_dict = build_local_dict(parameters)

    def _parse_string(expr_str: str):
        expr_str = expr_str.strip()
        return parse_expr(
            expr_str,
            transformations=transformations,
            local_dict=local_dict,
            evaluate=True,
        )

    try:
        lower = _parse_string(lower_str)
        upper = _parse_string(upper_str)
        return [lower, upper]

    except Exception as e:
        raise ValueError(f"Invalid bounds: {bounds}\nError: {e}")


def find_constants(parametrization, parameters):
    """
    Identify arbitrary constants in parametrization
    """
    if not isinstance(parameters, list):
        parameters = [parameters]

    param_set = {str(p) for p in parameters}

    quantity_set = set()
    symbol_set = set()
    constants = []

    for expr in parametrization:
        if not isinstance(expr, Expr):
            continue

        quantity_set |= expr.atoms(Quantity)
        symbol_set |= expr.free_symbols

    for symbol in symbol_set:
        if str(symbol) not in param_set:
            constants.append(Symbol(str(symbol), real=True, positive=True))

    return constants
