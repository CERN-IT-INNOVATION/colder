

import numpy as np
import sympy
import inspect

from typing import Tuple, List, Union

import colder.core.pauli_algebra as palg
import colder.simulation.numerical as cnum
from colder.core.subroutines import *





class hamiltonian:
    def __init__(self, operator : str, targets : list, coeff : str, 
            coeff_function : callable = None, target_coeffs : Union[np.ndarray,None,int,float] = None
        ) -> None:
        """
        Hamiltonian object.

        Args:
            operator (str): String defining the Pauli operators to apply.
            targets (list): List of spin targets to apply the operators. The length of the tuples has to match the length of `operator` string.
            coeff (str): Coefficient symbol.
            coeff_function (callable, optional): Function to associate to coefficient. To use COLD simulation, this argument has to be provided. Defaults to None.
            target_coeffs (Union[np.ndarray,None,int,float], optional): Coefficients to associate to the target spins. Defaults to None.

        Raises:
            Exception: First argument of `coeff_function`, if provided, must be `t` (time dependence).
        """
        self.targets : list = targets     # list of tuples
        if not isinstance(self.targets[0], tuple):
            self.targets = [ tuple([e]) for e in self.targets ]
        self.operator : str = operator    # string of operators to apply for each tuple in target site list
        
        self.n_terms : int = len(self.targets)          # number of targets
        self.operator_order : int = len(self.operator)  # how many operators are multiplied in a single tensor product
        
        # inject coefficient name and callable functions
        self.coeff : str = coeff
        self.coeff_function : callable = coeff_function
        
        if self.coeff_function is not None:
            arg = self.get_coeff_function_arguments()
            if arg[0] != 't':  # TODO: is this really necessary?
                raise Exception('first argument must be t')
        
        # process target coefficient
        if target_coeffs is not None:
            if hasattr(target_coeffs, '__iter__'):
                # if it is iterable, check for len
                assert len(target_coeffs) == self.n_terms
            else:
                # if it is not iterable, broadcast to np.ndarray
                target_coeffs = np.ones(self.n_terms)*target_coeffs 
            
        self.target_coeffs : Union[np.array,None] = target_coeffs

    
    
    def expression(self, sum_symbol : str = 'i', weight_symbol : str = 'w', operator_symbol : str = 'sigma_'):
        """
        Print a dummy expression for current hamiltonian.

        Args:
            sum_symbol (str, optional): Symbol to use in the sum. Defaults to `i`.
            weight_symbol (str, optional): Placeholder letter for symbols. Defaults to `w`.
            operator_symbol (str, optional): Spin operator symbol. Defaults to `sigma_`.
        """
        sum_expression = operator_symbol + self.operator + '^i'
        global_coeff = None
        
        # add weights, if any
        if self.target_coeffs is not None:
            if np.all(self.target_coeffs == self.target_coeffs[0]):
                # all elements are equal
                global_coeff : float = self.target_coeffs[0]
            else:
                sum_expression = weight_symbol + '_' + sum_symbol + '*' + sum_expression
        
        if global_coeff is not None:
            return sympy.Symbol(self.coeff) *global_coeff* sympy.Sum(sum_expression, (sympy.Symbol(sum_symbol), 1, len(self.targets)) )
        else:
            return sympy.Symbol(self.coeff) * sympy.Sum(sum_expression, (sympy.Symbol(sum_symbol), 1, len(self.targets)) )
    
    
    def get_strings(self, total_sites : int) -> dict[str, complex]:
        """
        Returns full strings with coefficients from current hamiltonian.

        Args:
            total_sites (int): Number of total spin sites.

        Returns:
            dict[str, complex]: Dictionary of expressions and coefficients.
        """
        if self.target_coeffs is None:   w = np.ones(self.n_terms)
        else:                            w = self.target_coeffs
        return { palg.build_string_from_operator_and_target(strlen=total_sites, target=tt, operator_sequence=self.operator) : cc for tt, cc in zip(self.targets, w) }
    
    
    def __repr__(self):
        return (
            f"{self.__class__.__name__}"
            f"(operator={self.operator}, coeff={self.coeff})"
        )
        
    
    def __add__(self, otherH):
        if otherH is None:  return self
        return hamiltonian_collection(self, otherH)
    
    def __iadd__(self, otherH):
        if otherH is None:  return self
        return hamiltonian_collection(self, otherH)
    
    
    def get_coeff(self) -> str:
        """
        Get the coefficient of this hamiltonian.

        Returns:
            str: Coefficient string.
        """
        return self.coeff
    
    
    def get_coeff_function_arguments(self, exception_if_none : bool = True) -> List[str]:
        """
        Get the argument name for the coefficient function.

        Args:
            exception_if_none (bool, optional): If True, raises an exception if there exist no function coefficient. Else, returns an empty list. Defaults to True.

        Raises:
            Exception: If `exception_if_none`, the exception is raised if no coefficient function is provided.

        Returns:
            List[str]: List of arguments for the coefficient function.
        """
        if self.coeff_function is None:
            if exception_if_none:
                raise Exception('coeff function not provided in __init__')
            else:
                return []
        return inspect.getfullargspec(self.coeff_function)[0]
    
    
    def make_coefficient_interpolation(self, trange : tuple, fargs : dict, n_interp_points : int = 50, **kwargs) -> cnum.interpolator1D:
        """
        Interpolates the coefficient function in time range `trange`.

        Args:
            trange (tuple): Range of time to interpolate.
            fargs (dict): Arguments to be passed to the coefficient function.
            n_interp_points (int, optional): Number of points to interpolate. Defaults to 50.

        Raises:
            Exception: Coefficient function must be provided to use this feature.

        Returns:
            cnum.interpolator1D: Numerical spline interpolator for coefficient function.
        """
        if self.coeff_function is None:
            raise Exception(f'coeff function {self.coeff} not provided in __init__')
        return cnum.interpolator1D.from_function(self.coeff_function, trange, n_interp_points = n_interp_points, fargs=fargs, **kwargs)






class hamiltonian_collection:
    """
    This class creates collection of hamiltonian object that have to be summed.

    Example:
        .. highlight:: python
        .. code-block:: python

        import colder.core.physics as cphys
        
        H_Z = cphys.hamiltonian('Z', [(0,) , (1,)], coeff = 'Z')
        H_X = cphys.hamiltonian('X', [(0,) , (1,)], coeff = 'X')      
        
        H = H_Z + H_X   # this is a collection (type: cphys.hamiltonian_collection)
    """
    
    def __init__(self, *args):
        """
        Creates an hamiltonian collection "summing" instances of hamiltonians.
        """
        self.terms = []
        
        for tt in locals()['args']:
            if isinstance(tt, hamiltonian_collection):
                # hamiltonian collections are unwrapped
                self.terms += tt.terms
            else:
                self.terms.append(tt)
    
    def get_coeffs(self) -> List[str]:
        return [ term.coeff for term in self.terms ]
    
    def get_unique_coeffs(self) -> List[str]:
        return list( set( self.get_coeffs() ) )
    
    def get_strings(self, total_sites : int) -> List[ dict[str, complex] ]:
        return [ term.get_strings(total_sites) for term in self.terms ]
    
    def expression(self, group_coefficients : bool = True):
        """
        Print a dummy expression for current hamiltonian.

        Args:
            group_coefficients (bool, optional): If True, coefficients are grouped, if needed. Defaults to True.
        """
        expr = sum( ee.expression() for ee in self.terms )
        if group_coefficients:
            expr = sympy.collect(expr, [ sympy.Symbol(ee) for ee in self.get_unique_coeffs() ] )
        return expr

    def make_symbolic_expression(self, total_sites : int):
        """
        Returns the symbolic expression for this hamiltonian collection.

        Args:
            total_sites (int): Total number of spin sites.
        """
        return sum( sympy.Symbol(cc)*make_sum_expression(st) for cc, st in zip(self.get_coeffs(), self.get_strings(total_sites)) )
    
    def __add__(self, otherH):
        return hamiltonian_collection(self, otherH)
    
    def __repr__(self):
        return (
            f"{self.__class__.__name__}"
            f"(terms={self.terms}, coeff={self.get_unique_coeffs()})"
        )
        
    def make_coefficients_interpolation(self, trange : tuple, fargs : dict, make_unique : bool = True, **kwargs) -> dict[str, cnum.interpolator1D]:
        """
        Interpolates the coefficient function in time range `trange` for each unique coefficient of hamiltonian collection.

        Args:
            trange (tuple): Range of time to interpolate.
            fargs (dict): Arguments to be passed to the coefficient function.
            make_unique (bool, optional): If True, uses unique coefficients to avoid multiple computation of the same interpolation. Defaults to True.

        Returns:
            dict[str, cnum.interpolator1D]: Dictionary of interpolated functions.
        """
        if make_unique:
            # strongly suggested, to avoid interpolating many times the same parameter
            unique_coeff_objects : List[str] = { hh.coeff : hh for hh in self.terms }
            return { kk : hh.make_coefficient_interpolation(trange=trange, fargs=fargs, **kwargs) for kk, hh in unique_coeff_objects.items() }
        
        return { hh.coeff : hh.make_coefficient_interpolation(trange=trange, fargs=fargs, **kwargs) for hh in self.terms }





class empty_hamiltonian:
    """
    Empty hamiltonian class. The purpuse of this dummy class is to provide a custom empty object to be summed with hamiltonians.
    
    Example:
        .. highlight:: python
        .. code-block:: python

            import colder.core.physics as cphys
            
            H = cphys.empty_hamiltonian()
            
            # just a terrible example:
            for i in range(3):
                H += cphys.hamiltonian('XX', [(i,i+1)], coeff = 'J', coeff_function=Jf, target_coeffs= 1/(i+1) )

    """
    def __init__(self) -> None:
        pass
    
    def __repr__(self):
        return (
            f"{self.__class__.__name__}"
            f" use + or += operator to add a hamiltonian obj"
        )
    
    def __add__(self, other):
        return other
    
    def __iadd__(self, other):
        return other
