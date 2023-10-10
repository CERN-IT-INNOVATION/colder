
import numpy as np
import scipy

import qibo

from typing import Tuple, List, Union


def get_lower_eigen(hh, k : int = 8) -> Union[np.ndarray, np.ndarray]:
    assert k>0, 'must take at least one eigenvalues'
    return evals, evect


def get_superposition(hh, n_eig : int) -> np.ndarray:
    raise Exception('TODO, not yet implemented')
    return psi


def get_groundstate_superposition(hh, n_eig : int = 8, method : str = 'dense', print_info : bool = False) -> np.ndarray:
    
    assert method in ['sym', 'dense'], 'unknown eigenvalue routine method'
    
    if method == 'sym':
        # ISSUE -----------------------
        #   Qibo 0.1.15 calculates the dense form of a symbolic Hamiltonian anyway!
        #   So this way of computing the evals and evectors is not better than 'dense' method...
        evals, evect = hh.eigenvalues(k=n_eig), hh.eigenvectors(k=n_eig)
    
    elif method == 'dense':
        matrix = hh.matrix 
        evals, evect = scipy.linalg.eigh(matrix, subset_by_index = (0,n_eig))

    mask : np.ndarray = (evals == evals[0])
    n_superposition : int = np.sum(mask)
    gs : np.ndarray = np.sum( evect[:,mask].T, axis = 0)/np.sqrt(n_superposition)
    
    if n_superposition == n_eig:
        print('warning: Number of selected eigenvectors is equal to max allowed eigenvalues. There might be more states with this eigenvalue.')
    if print_info:
        print('[info] gs is superposition of {} states'.format(n_superposition) )
    
    return gs



def timedependent_evolution(htf : callable, T : float, P : int, 
        psi0 : Union[np.ndarray, str] = 'superpos', gs_args : dict = {},
        hook_function : None = None
    ) -> np.ndarray:
    """Execute a time evolution via Qibo state evolution. The Hamiltonian is a callable function of time H(t)."""
    
    if hook_function is not None:
        raise Exception('hook function is not supported with qibo backend')
    
    # compute ground state of ht(0) to initialize the state, if not provided
    if isinstance(psi0, str):
        if psi0 == 'superpos':
            psi0 = get_groundstate_superposition( htf(0) )
        elif psi0 == 'qibo':
            psi0 = htf(0).ground_state()
        else:
            raise Exception(psi0 + ' is not valid string identifier for psi0')
    
    # evolve
    evolve = qibo.models.StateEvolution(htf, dt = T/P)
    final_state = evolve(final_time=T, initial_state=psi0)
    
    return final_state