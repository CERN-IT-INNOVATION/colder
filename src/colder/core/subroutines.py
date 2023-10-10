
import sympy

from typing import Tuple, List



# GENERIC DATA MANIPULATION

def remove_null_values_dictionary(expr : dict) -> dict:
    """
    Removes entries with null coefficients from expression dictionary.

    Args:
        expr (dict): Dictionary of strings and related coefficients.

    Returns:
        dict: Dictionary without null values.
    """
    return { x : y for x, y in expr.items() if y!=0 }




# STRING BUILDERS

def build_string_from_operator_and_target(strlen : int, target : tuple, operator_sequence : str, default_op : str = 'I') -> str:
    """
    Create a string inserting the operator_sequence operators in target indices.
    
    For instance, calling the function with arguments  5, (1,3), 'XY' will return 'IXIYII'.

    Args:
        strlen (int): Total length of the string
        target (tuple): Indices of string to replace
        operator_sequence (str): Operators to replace (shoul have the same number of characters as target tuple elements)
        default_op (str, optional): default single site operator. Defaults to 'I' (identity in the algebra).

    Returns:
        str: String with operator in target position.
    """
    this_string : list = [default_op] * strlen
    for ii, op in zip(target, operator_sequence):  this_string[ii] = op
    return ''.join(this_string)


def build_string_from_regular_pattern(N : int, pattern : str, default_op : str = 'I', coeff : complex = 1) -> dict[str, complex]:
    """
    Given a string pattern, builds all the strings in which the pattern is repeated through all the N sites.

    Args:
        N (int): maximum length of string
        pattern (str): pattern (of length L) to be repeated
        default_op (str, optional): Default operator to fill the string. Defaults to 'I', for identity operator.
        coeff (complex, optional): Default coefficient to assign. Defaults to 1.

    Returns:
        dict[str, complex]: Dictionary of strings and related coefficients.
    """
    L = len(pattern)
    assert N >= L, 'pattern len is larger than system size'
    
    return { default_op*i + pattern + default_op*(N-L-i) : coeff for i in range(N + 1 - L) }





# EXPRESSION CREATION

def make_sum_expression(ops : dict[str, complex], global_coeff : complex = 1, symbol_prefix : str = 'sigma_'):
    """Sum the operators in input dictionary with related coefficients.

    Args:
        ops (dict[str, complex]): Strings to be cast as operators
        coeff (complex, optional): Global coefficient. Defaults to 1.
        
    Returns:
        _type_: _description_
    """
    # TODO
    return sum( global_coeff*y*sympy.Symbol(symbol_prefix+x) if global_coeff*y != 1 else sympy.Symbol(symbol_prefix+x) for x,y in ops.items() )




# EXPRESSION MANIPULATION


def expression_singular_square(expr, collectors):
    tt = sympy.collect(expr, collectors, evaluate = False)
    
    new_expr = None
    for op, ee in tt.items():
        if new_expr is None:
            new_expr = op*(ee**2)
        else:
            new_expr += op*(ee**2)
    return new_expr


def expressions_to_linear_system_matrix( expressions : list, variables : list) -> sympy.matrices:
    
    # create matrix A
    A = sympy.Matrix( [[ee.coeff( sympy.Symbol(var) ) for var in variables] for ee in expressions] )
    # create matrix B
    B = sympy.Matrix([-ee + sum(ee.coeff( sympy.Symbol(var) ) * sympy.Symbol(var) for var in variables) for ee in expressions ])
    return A, B

def make_linear_system(A : sympy.matrices, B : sympy.matrices, variables) -> sympy.Eq:
    return sympy.Eq(A * sympy.Matrix(variables), B)



def get_operators_from_expr_with_prefix(expr, operator_prefix : str) -> List:
    """Returns a list of operators (as Sympy symbols) in expr starting with operator_prefix."""
    symbols = expr.free_symbols
    return [ element for element in symbols if element.name.startswith(operator_prefix) ]

def replace_operators_in_expr(expr, operator_prefix : str, replacement_value : float = 1):
    sigmas = get_operators_from_expr_with_prefix(expr, operator_prefix)
    for ss in sigmas: expr = expr.subs(ss, replacement_value)
    return expr



def get_free_symbols_name(expr) -> set:
    syms = [ ss.name for ss in expr.free_symbols ]
    return set( syms )