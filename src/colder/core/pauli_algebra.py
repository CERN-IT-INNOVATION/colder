
from colder.core.subroutines import *

import itertools
from typing import Tuple, List




def pauli_product(a : int, b : int) -> Tuple[complex, int]:
    """
    Simplify the product between the Pauli operators sigma_a and sigma_b and compute the coefficient.
    The operators are ancoded using integer values:  1 -> sigma_x, 2 -> sigma_y, 3 -> sigma_z.
    0 is the Identity.

    Args:
        a (int): first Pauli operator
        b (int): second Pauli operator (order matters)

    Returns:
        Tuple[complex, int]: coefficient and resulting operator
    """
    
    # NOTE   operators are encoded as integers:
    #     0 ->  Identity
    #     1 ->  sigma_x    2 -> sigma_y    3 -> sigma_z
    
    # sigma_i * sigma_i = I
    if a == b:       return 1, 0
    
    # one of two terms is I
    if a*b == 0:     return 1, max(a,b)
    
    # nontrivial cases
    else:
        return 1j if (a,b) in [(1,2), (2,3), (3,1)] else -1j, 6 - a - b




def pauli_sites_product(A : str, B : str, coeff : complex = None) -> Tuple[complex, str]:
    """
    Given two input strings encoding Pauli operators in local sites, computes the local product of Pauli operators and returns the coefficient and resulting operator string.

    Args:
        A (str): String encoding Pauli operator (example: 'IIXY')
        B (str): _description_
        coeff (complex, optional): Global coefficient to be multiplied to the final result coefficient. Defaults to None.

    Returns:
        Tuple[complex, str]: The coefficient and the final operator string.
    """
    
    # checks: strings should be uppercase
    A = A.upper()
    B = B.upper()
    assert len(A) == len(B), 'input strings must have same number of sites'
    
    if coeff is None:  coeff = 1 + 0j
    
    # NOTE: Pauli operators are encoded as int using the absolute position on this string
    paulis : str = 'IXYZ'
    
    word : str = ''
    for a, b in zip(A, B):
        a, b = paulis.find(a), paulis.find(b)
        cc, op = pauli_product(a,b)
        coeff *= cc
        word += paulis[op]
        
    return coeff, word 



# from now on the term operator refers to a dictionary of strings (keys) and coefficients (value)


def operator_product(a : dict[str, complex], b : dict[str, complex], expr : dict = None, coeff : complex = +1, remove_null : bool = True) -> dict[str, complex]:
    """ 
    Compute product between string prompted as dictionaries.

    Args:
        a (dict[str, complex]): First set of strings.
        b (dict[str, complex]): Second set of strings (order matters).
        expr (dict, optional): Dictionary to expand with terms. Defaults to None.
        coeff (complex, optional): Global coefficient for strings. Defaults to +1.
        remove_null (bool, optional): Remove null terms. Defaults to True.
        
    Returns:
        dict[str, complex]: Dictionary of string and associated coefficients.
    """
    if expr is None: expr = {}
    # bug warning: do not use expr = {} in function declaration
    
    for x, y in itertools.product(a.items(), b.items()):
        # NOTE:  x[0] is the string, x[1] is the coefficient
        #        and likewise for y
        c, s = pauli_sites_product(x[0], y[0])
        c *= x[1]*y[1]
        
        if s not in expr:
            expr[s] = coeff*c
        else:
            expr[s] += coeff*c
    
    if remove_null:  expr = remove_null_values_dictionary(expr)
    return expr






def operator_commutator(a : dict[str, complex], b : dict[str, complex], coeff : complex = 1, remove_null : bool = True) -> dict[str, complex]:
    """
    Compute commutator (ab-ba) between strings prompted as dictionaries.


    Args:
        a (dict[str, complex]): _description_
        b (dict[str, complex]): _description_
        coeff (complex, optional): _description_. Defaults to 1.
        remove_null (bool, optional): _description_. Defaults to True.

    Returns:
        dict[str, complex]: _description_
    """
    expr : dict = operator_product(a, b, coeff=+1, remove_null=False)
    # not removing here to facilitate search in next term
    expr = operator_product(b, a, expr=expr, coeff=-1, remove_null=remove_null)
    
    # multiply global coefficient
    if coeff != 1:   expr = { x : coeff*y for x,y in expr.items() }
    
    if remove_null:  expr = remove_null_values_dictionary(expr)
    return expr

