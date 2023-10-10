
import numpy as np
import math

import qibo
import qibo.symbols as qsy
from colder.core.physics import hamiltonian, hamiltonian_collection

from typing import Tuple, List, Union


# interface with physics.hamiltonian
class interface:
    
    __operator_associations = { 'X' : qsy.X, 'Y' : qsy.Y, 'Z' : qsy.Z }
    
    def __init__(self, nspin : int, hh : Union[hamiltonian, hamiltonian_collection]):
        
        # store internally the target hamiltonian objects
        self.hobjs = []
        self.nspin = nspin 
        
        if isinstance(hh, hamiltonian_collection):
            # hamiltonian collections are unwrapped
            self.hobjs += hh.terms
        elif isinstance(hh, hamiltonian):
            # hamiltonian terms are just appended
            self.hobjs.append(hh)
        else:
            raise Exception('unknown input type')
        
    
    def __evaluate_operator_term(self, hh : hamiltonian, idx : int):
        """Compute the qibo symbolic operator for term idx of hamiltonian hh, i.e. terms"""
        
        # retrieve the qibo operators associated to this hamiltonian
        opf : list = [ self.__operator_associations[ii] for ii in hh.operator ]
        
        # get the current term spin indices
        sites : tuple = tuple( hh.targets[idx] )
        assert len(sites) == len(opf), 'hamiltonian operators and site tuple must have same length'
        
        return math.prod( fun(site) for fun, site in zip(opf, sites) )
    
    def __evaluate_static_terms(self, hh : hamiltonian):
        """Computes the static term for a (single) hamiltonian object. Sum(coeff*terms)"""
        
        # retrieve the coefficients (and patch for ones)
        target_coeffs = hh.target_coeffs
        if hh.target_coeffs is None:  target_coeffs = np.ones(hh.n_terms)
        
        return sum( target_coeffs[ii]*self.__evaluate_operator_term(hh, ii) for ii in range(hh.n_terms) )
    
        
    def compute_timeindependent_expressions(self) -> List:
        return [ self.__evaluate_static_terms(thishh) for thishh in self.hobjs ]
    
    def compute_timeindependent_expressions_with_coeff(self) -> List:
        return [ (hh.coeff, self.__evaluate_static_terms(hh)) for hh in self.hobjs ]
    
    def evaluate_static_qibo_hamiltonians(self, nqubits : int, dense : bool = False) -> List: # TODO deprecate
        raise Exception('should not be used')
        return [ qibo.hamiltonians.SymbolicHamiltonian(ee, nqubits = nqubits, dense=dense) for ee in self.evaluate_static_qibo_expressions() ]
    
    def finalize(self, terms):
        # action to apply before making it a function
        return qibo.hamiltonians.SymbolicHamiltonian(terms, nqubits = self.nspin)

    
    # this is required by COLD interface ----------
    def COLD_make_hamiltonian(self):
        return self.compute_timeindependent_expressions_with_coeff()

    def COLD_get_finalize_method(self):
        def finalize(tt):
            return qibo.hamiltonians.SymbolicHamiltonian(tt, nqubits = self.nspin)
        return finalize