from sympy import *
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication_application,
    convert_xor,
    function_exponentiation,
)

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

allowed_symbols = {
    name: symbols(name, real=True) for name in (english_letters + greek_letters)
}

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

local_dict = {}
local_dict.update(allowed_symbols)
local_dict.update(allowed_functions)
local_dict.update(allowed_constants)

transformations = (
    standard_transformations
    + (implicit_multiplication_application,)
    + (convert_xor,)
    + (function_exponentiation,)
)


def parse_input(expr_string: str):
    """
    Parse plaintext input into a sympy expression
    """

    # Remove whitespace
    expr_string = expr_string.strip()

    try:
        parsed = parse_expr(
            expr_string,
            transformations="all",
            local_dict=local_dict,
            evaluate=True,
        )

        return parsed
    except Exception as e:
        raise ValueError(f"Invalid expression: {expr_string}\nError: {e}")
