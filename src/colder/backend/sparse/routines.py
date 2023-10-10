
import numpy as np
import scipy

from typing import Tuple, List, Union


def get_lower_eigen(hh, k : int = 8, method : str = 'dense'):
    assert k>0, 'must take at least one eigenvalues'
    assert method in ['sparse', 'dense'], 'unknown eigenvalue routine method'
    
    if method == 'sparse':
        # ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html
        evals, evect = scipy.sparse.linalg.eigsh(hh, k=k, which='SA', tol = 10**-8)
        
    elif method == 'dense':
        matrix = hh.todense()
        evals, evect = scipy.linalg.eigh(hh, subset_by_index = (0,k-1))
        
    return evals, evect


def get_superposition(hh, n_eig : int, method : str = 'dense'):
    
    if method == 'sparse':
        # ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html
        evals, evect = scipy.sparse.linalg.eigsh(hh, k=n_eig-1, which='SA', tol = 10**-8)
        
    elif method == 'dense':
        matrix = hh.todense()
        evals, evect = scipy.linalg.eigh(hh, subset_by_index = (0,n_eig-1))
        
    psi : np.ndarray = np.sum( evect[:,:n_eig].T, axis = 0)/np.sqrt(n_eig)
    return psi



def get_groundstate_superposition(hh, n_eig : int = 8, method : str = 'dense', print_info : bool = False):
    
    assert method in ['sparse', 'dense'], 'unknown eigenvalue routine method'
    
    if method == 'sparse':
        # ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html
        evals, evect = scipy.sparse.linalg.eigsh(hh, k=n_eig-1, which='SA', tol = 10**-8)
        
    elif method == 'dense':
        matrix = hh.todense()
        evals, evect = scipy.linalg.eigh(matrix, subset_by_index = (0,n_eig-1))

    mask : np.ndarray = (evals == evals[0])
    n_superposition : int = np.sum(mask)
    gs : np.ndarray = np.sum( evect[:,mask].T, axis = 0)/np.sqrt(n_superposition)
    
    if (n_superposition == n_eig) and (n_eig != 1):
        print('warning: Number of selected eigenvectors is equal to max allowed eigenvalues. There might be more states with this eigenvalue.')
    if print_info:
        print('[info] gs is superposition of {} states'.format(n_superposition) )
    
    return gs



def timedependent_evolution(htf : callable, T : float, P : int, psi0 : Union[np.ndarray, str] = 'superpos', gs_args : dict = {}, force_norm : bool = True, hook_function : Union[callable,None] = None):
    """Execute a time evolution. The Hamiltonian is a callable function of time H(t)."""
    
    # compute ground state of ht(0) to initialize the state, if not provided
    if isinstance(psi0, str):
        if psi0 == 'superpos':
            psi0 = get_groundstate_superposition( htf(0), *gs_args )
        elif psi0 == 'single':
            psi0 = get_groundstate_superposition( htf(0), n_eig=1, *gs_args )
        else:
            raise Exception(psi0 + ' is not valid string identifier for psi0')
    
    state = psi0
    time = np.linspace(0,T,P+1)
    dt = np.mean( np.diff(time) )
    
    # FIXME
    
    def _single_step_expm(Ht, state):
        U = scipy.sparse.linalg.expm(-1.j*dt*Ht)
        return U * state
    
    def _single_step_sparse(Ht, state):
        evals, evect = scipy.sparse.linalg.eigsh(Ht, k = 30, which = 'SA')
        U = np.exp( -1.j*evals*dt )
        return np.dot(evect, U * np.dot(evect.transpose().conjugate(), state))
    
    def _single_step_dense(Ht, state):
        evals, evect = scipy.linalg.eigh( Ht.todense() )
        U = np.exp( -1.j*evals*dt )
        return np.dot(evect, U * np.dot(evect.transpose().conjugate(), state))

    if hook_function is None:
        # loop without hook function
        for t in time[1:]:
            state = _single_step_dense( htf(t), state)
            if force_norm: state = state/(state @ state.conj())
            
    else:
        # loop with hook function call, exactly the same as above for the rest
        for t in time[1:]:
            state = _single_step_expm( htf(t), state)
            if force_norm: state = state/(state @ state.conj())
            hook_function(t, state)
            
    return state