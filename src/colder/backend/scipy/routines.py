
import numpy as np
import scipy

from typing import Tuple, List, Union




def get_lower_eigen(hh, k : int = 8) -> Union[np.ndarray, np.ndarray]:
    assert k>0, 'must take at least one eigenvalues'
    
    evals, evect = scipy.linalg.eigh(hh, subset_by_index = (0,k-1))
    return evals, evect


def get_superposition(hh, n_eig : int) -> np.ndarray:
    
    evals, evect = scipy.linalg.eigh(hh, subset_by_index = (0,n_eig))
    psi : np.ndarray = np.sum( evect[:,:n_eig].T, axis = 0)/np.sqrt(n_eig)
    return psi


def get_groundstate_superposition(hh, n_eig : int = 8, print_info : bool = False) -> np.ndarray:
    
    evals, evect = scipy.linalg.eigh(hh, subset_by_index = (0,n_eig-1))

    mask : np.ndarray = (evals == evals[0])
    n_superposition : int = np.sum(mask)
    gs : np.ndarray = np.sum( evect[:,mask].T, axis = 0)/np.sqrt(n_superposition)
    
    if (n_superposition == n_eig) and (n_eig != 1):
        print('warning: Number of selected eigenvectors is equal to max allowed eigenvalues. There might be more states with this eigenvalue.')
    if print_info:
        print('[info] gs is superposition of {} states'.format(n_superposition) )
    
    return gs



def timedependent_evolution(htf : callable, T : float, P : int, 
    psi0 : Union[np.ndarray, str] = 'superpos', gs_args : dict = {},
    force_norm : bool = True, hook_function : Union[callable,None] = None
    ) -> np.ndarray:
    """Execute a time evolution. The Hamiltonian is a callable function of time H(t)."""
    
    # compute ground state of ht(0) to initialize the state, if not provided
    if isinstance(psi0, str):
        if psi0 == 'superpos':
            psi0 = get_groundstate_superposition( htf(0), *gs_args )
        elif psi0 == 'single':
            psi0 = get_groundstate_superposition( htf(0), n_eig=1, *gs_args )
        else:
            raise Exception(psi0 + ' is not valid string identifier for psi0')
    
    state : np.ndarray = psi0
    time = np.linspace(0,T,P+1)
    dt = np.mean( np.diff(time) )
    
    # these two methods should be equivalent
    def _single_step_expm(Ht, state):
        U = scipy.linalg.expm(-1.j*dt*Ht)
        return U @ state
    
    def _single_step_dense(Ht, state):
        evals, evect = scipy.linalg.eigh( Ht )
        U = np.exp( -1.j*evals*dt )
        return np.dot(evect, U * np.dot(evect.transpose().conjugate(), state))
    
    if hook_function is None:
        # loop without hook function
        for t in time[1:]:
            state = _single_step_expm( htf(t), state)
            # force normalization
            if force_norm: state = state/(state @ state.conj())
            
    else:
        # loop with hook function call, exactly the same as above for the rest
        for t in time[1:]:
            state = _single_step_expm( htf(t), state)
            if force_norm: state = state/(state @ state.conj())
            hook_function(t = t, state = state)
            
    return state