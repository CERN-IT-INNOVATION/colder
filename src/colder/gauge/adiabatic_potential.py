
import numpy as np
import matplotlib.pyplot as plt

import sympy
import itertools

import colder.core.physics as cphys
import colder.core.subroutines as csub
import colder.core.pauli_algebra as cpal
import colder.simulation.numerical as cnum


from typing import Tuple, List, Union




class cache_lcd_numpy:
    """
    This class collects the intermidiate result of AGP to be cached.
    """
    def __init__(self, agpobj):
        """
        Makes the LCD cache from AGP.
        
        Args:
            agpobj (adiabaticgp): adiabaticgp object.
            
        """
        A, B, syms = agpobj.make_lcd_numpy(modules = 'numpy')
        self.A = A
        self.B = B
        self.syms = syms


class adiabaticgp:
    def __init__(self,
            total_sites : int, H : cphys.hamiltonian_collection, ansatz : cphys.hamiltonian_collection,
            sort_params : bool = True, 
            enable_cache : bool = True, inject_cache : Union[None,cache_lcd_numpy] = None
        ) -> None:
        """
        Creates an Adiabatic Gauge Potential for system `H` and GP `ansatz`. The total number of spins must be provided.

        Args:
            total_sites (int): Number of sites.
            H (colder.core.physics.hamiltonian_collection): System hamiltonian.
            ansatz (colder.core.physics.hamiltonian_collection): Ansatz hamiltonian.
            sort_params (bool, optional): Sort parameters in alphabetic order. Defaults to True.
            enable_cache (bool, optional): If True, caching is enabled. Defaults to True.
            inject_cache (Union[None,cache_lcd_numpy], optional): If not None, use a cache object instead of computing the AGP from scratch. Defaults to None.

        Raises:
            Exception: The APG cache object is not valid.
        """
        self.total_sites = total_sites
        
        self.Hobj = H
        self.ansatzobj = ansatz
        
        # get unique parameters from ansatz
        self.model_unique_parameters = self.Hobj.get_unique_coeffs()
        self.ansatz_unique_parameters = self.ansatzobj.get_unique_coeffs()
        
        if sort_params:
            self.model_unique_parameters.sort()
            self.ansatz_unique_parameters.sort()
        
        self.operators_prefix = 'sigma_'
        self.derivative_prefix = 'd'
        
        # store last lambdify 
        self.__last_lambdify_args = None
        
        # cache for agp functions
        self.cache = enable_cache
        self.cache_lcd = None
        if enable_cache:
            if inject_cache is None: 
                self.make_cache()
            elif isinstance(inject_cache, cache_lcd_numpy):
                self.cache_lcd = inject_cache
            else:
                raise Exception('not valid agp cache object')
    
    
    def make_cache(self) -> None:
        """
        Computing the AGP and store it in cache.
        """
        print('[info] caching agp')
        self.cache_lcd = cache_lcd_numpy(self)
        self.cache = True
    
    def clear_cache(self) -> None:
        """Remove cached LCD. On next call, the AGP will be computed from scratch."""
        self.cache_lcd = None
        
    def set_cache(self, cc : cache_lcd_numpy) -> None:
        """
        Inject a cache object.

        Args:
            cc (cache_lcd_numpy): AGP cache.
        """
        self.cache = True
        self.cache_lcd = cc
    
    def get_cache(self) -> Union[cache_lcd_numpy, None]:
        """Returns the AGP cache, if computed.

        Returns:
            Union[cache_lcd_numpy, None]: AGP cache or None if not yet computed.
        """
        return self.cache_lcd


    def compute(self):
        """
        Compute simbolically :math:`G` as expression.
        """
        H = self.Hobj.make_symbolic_expression(self.total_sites)
        
        # compute the derivative by replacing the symbols with the dotted ones
        G = - H.subs( { sympy.Symbol(ss) : sympy.Symbol(self.derivative_prefix + ss) for ss in self.model_unique_parameters } )
        
        group_model = zip( self.Hobj.get_strings(self.total_sites), self.Hobj.get_coeffs() )
        group_ansatz = zip( self.ansatzobj.get_strings(self.total_sites), self.ansatzobj.get_coeffs() )
        
        for mm, aa in itertools.product(group_model, group_ansatz):
            comms = cpal.operator_commutator(mm[0], aa[0], coeff = 1j)
            G += sympy.Symbol(mm[1])*sympy.Symbol(aa[1])*csub.make_sum_expression(comms)
        
        return G.expand()
    
    def compute_square_trace(self):
        """
        Compute simbolically :math:`Tr[G^2]`.
        """
        G = self.compute()
        sigmas = csub.get_operators_from_expr_with_prefix(G, self.operators_prefix)
        trG2 = csub.expression_singular_square(G, sigmas)
        return csub.replace_operators_in_expr(trG2, self.operators_prefix, replacement_value = 1)
    
    
    def get_lcd_equation_byvar(self, var : str):
        """Computes the LCD equations for variable `var`.

        Args:
            var (str): Name of the variable to be collected.
        """
        trG2 = self.compute_square_trace()
        alpha_diff = sympy.diff( trG2, sympy.Symbol(var) ).expand()
        return sympy.collect(alpha_diff, [ sympy.Symbol(ee) for ee in self.ansatz_unique_parameters ])
    
    def get_lcd_equations(self) -> List:
        """
        Compute the LCD equations for all variables.

        Returns:
            List: List of symbolic equations for each variable in ansatz.
        """
        # NOTE: the order of variables is given by obj.ansatz_unique_parameters
        return [ self.get_lcd_equation_byvar(vv) for vv in self.ansatz_unique_parameters ]
    
    
    def make_lcd_system(self):
        """
        Compute simbolically the LCD equations as a linear system.
        """
        # compute the equations
        eqs : list = self.get_lcd_equations()
        
        # get the linear system as sympy matrices
        A, B = csub.expressions_to_linear_system_matrix(eqs, self.ansatz_unique_parameters)
        return A, B
    
    
    def retrieve_model_symbols(self) -> List:
        """
        Returns the list of all symbols for this AGP model.

        Returns:
            List: List of symbols in this AGP model.
        """
        expected_symbols : list = self.model_unique_parameters + [ self.derivative_prefix + ee for ee in self.model_unique_parameters ]
        return expected_symbols

    
    def make_lcd_numpy(self, symbols_check : bool = True, modules : str = 'numpy', args_override : List = None
        ) -> Union[callable, callable, List]:
        """
        Compute simbolically the LCD equations as a linear system and lambdify the matrices to numpy array depending on system symbols.
        Remark: this is the cached object.
        
        Args:
            symbols_check (bool, optional): If True, symbol consistency is checked. Defaults to True.
            modules (str, optional): Module for lambidify target. Defaults to 'numpy'. Strongly suggest to not change this, unless you know what you are doing.
            args_override (List, optional): Override `expected_symbols`. Defaults to None.

        Raises:
            Exception: Unexpected symbols in equation.

        Returns:
            Union[callable, callable, List]: Callable functions for A, B. List of symbols to provide as arguments to A and B.
        """
        # get the linear system as matrices
        A, B = self.make_lcd_system()
        
        if args_override is None:
            # get (all) expected symbols from model parameter
            expected_symbols : list = self.retrieve_model_symbols()
        else:
            # override list of symbols
            expected_symbols = args_override
        
        if symbols_check:
            # get all the symbols from linear system
            symbols_eqs : set = csub.get_free_symbols_name( A )
            symbols_eqs.update( csub.get_free_symbols_name(B) )
            
            # check that all expected symbols are in equation symbols
            check = symbols_eqs.issubset( set(expected_symbols) )
            if not check:
                raise Exception('expression has unexpected symbols')

        # lambdify the matrices
        A_func = sympy.lambdify(expected_symbols, A, modules=modules)
        B_func = sympy.lambdify(expected_symbols, B, modules=modules)
        
        # store the last lambdify arguments, for redundancy
        self.__last_lambdify_args = expected_symbols
        
        return A_func, B_func, expected_symbols
        # note: the result of this function is cached
    
    
    
    def lcd_numerical_solver(self, time : np.ndarray, finj : dict, fargs : dict = {}, nocache : bool = False) -> np.ndarray:
        """
        Solve LCD equations for `time`, given a schedule injected through `finj`.

        Args:
            time (np.ndarray): Time array to solve for.
            finj (dict): Dictionary of callable functions to be called for each time step.
            fargs (dict, optional): Arguments of the interpolator functions. Defaults to {}.
            nocache (bool, optional): If True, forces the computation of LCD from scratch. Otherwise, cache is used, if available. Defaults to False.

        Raises:
            Exception: Injection dictionary has unexpected arguments. Are you sure the cache is correct?

        Returns:
            np.ndarray: Numerical solution for each unique parameter of this model.
        """
        if nocache:
            A, B, syms = self.make_lcd_numpy(modules = 'numpy')
        else:
            if self.cache:
                if self.cache_lcd is None:
                    self.make_cache()
                
                cc = self.cache_lcd
                A, B, syms = cc.A, cc.B, cc.syms
            else:
                A, B, syms = self.make_lcd_numpy(modules = 'numpy')
            
        # check that finj has all therequired arguments
        if set(syms) != set( finj.keys() ):
            raise Exception('injection dictionary has not the same keys expected by A B functions')
        
        numsol : np.matrix = np.empty( (len(time),len(self.ansatz_unique_parameters)) )
        for ii, tt in enumerate(time):
            
            time_args = { k : v(tt, **fargs) for k, v in finj.items() } 
            A_timed = A( **time_args )
            B_timed = B( **time_args )
            
            numsol[ii,:] = np.linalg.solve(A_timed, B_timed).flatten()
        
        return numsol
    
    
    def lcd_numerical_solver_from_range(self, trange : tuple, nsamples : int, finj : dict, fargs : dict = {}) -> np.ndarray:
        """
        Call `lcd_numerical_solver` for a time array in range `trange`.

        Args:
            trange (tuple): Range of time to solve for.
            nsamples (int): Number of time steps.
            finj (dict): Dictionary of callable functions to be called for each time step.
            fargs (dict, optional): Arguments of the interpolator functions. Defaults to {}.

        Returns:
            np.ndarray: Numerical solution for each unique parameter of this model.
        """
        time : np.array = np.linspace(trange[0], trange[1], nsamples)
        return self.lcd_numerical_solver(time, finj=finj, fargs=fargs)
    
    
    def lcd_interpolated_solver(self, trange : tuple, nsamples : int, finj : dict, fargs : dict = {}, bc_type : str = 'clamped') -> dict[str, cnum.interpolator1D]:
        """
        Solve numerically the LCD for time in `trange` and return a dictionary of values for each unique parameter of the system.

        Args:
            trange (tuple): Range of time to solve.
            nsamples (int): Number of time samples.
            finj (dict): Dictionary of callable functions to be called for each time step.
            fargs (dict, optional): Arguments of the interpolator functions. Defaults to {}.
            bc_type (str, optional): Spline option. Defaults to 'clamped'.

        Returns:
            dict[str, cnum.interpolator1D]: Dictionary for each interpolator associated to unique driving parameters of the system.
        """
        time : np.array = np.linspace(trange[0], trange[1], nsamples)
        sol : np.matrix = self.lcd_numerical_solver(time, finj=finj, fargs=fargs)
        
        solvers : dict = { k : None for k in self.ansatz_unique_parameters }
        
        # loop over every variable and interpolate
        for variable, yyy in zip(self.ansatz_unique_parameters, sol.T):
            solvers[variable] = cnum.interpolator1D(time, yyy, bc_type=bc_type)
        
        return solvers
    
    
    def match_interpolated_injections(self, system_interpolators : dict) -> dict:
        """
        Returns a dictionary for of callable function for the parameters in system and their derivatives.

        Args:
            system_interpolators (dict): Dictionary of interpolators.

        Returns:
            dict: Dictionary of interpolators for function and derivatives.
        """
        finj : dict = {}
        for coeff, interp in system_interpolators.items():
            finj[coeff] = interp.get_function()
            finj[self.derivative_prefix+coeff] = interp.get_derivative()
        
        return finj







def plot_interpolators(data : dict, time : np.array):
    """
    Plot the interpolated data in time array.

    Args:
        data (dict): Dictionary of arrays.
        time (np.array): Time values.
    """
    for k, interp in data.items():
        plt.plot(time, interp.get_function()(time), label = k)
