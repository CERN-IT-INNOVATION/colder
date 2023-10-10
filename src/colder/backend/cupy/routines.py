
import cupy
import cupyx

from typing import Tuple, List, Union



def get_lower_eigen(hh, k : int = 8) -> Union[cupy.ndarray, cupy.ndarray]:
    assert k>0, 'must take at least one eigenvalues'
    
    evals, evect = cupy.linalg.eigh(hh)
    return evals[:k].get(), evect[:,:k].get()


def get_superposition(hh, n_eig : int) -> cupy.ndarray:
    
    evals, evect = cupy.linalg.eigh(hh)
    psi : cupy.ndarray = cupy.sum( evect[:,:n_eig].T, axis = 0)/cupy.sqrt(n_eig)
    return psi



def get_groundstate_superposition(hh, n_eig : int = 8, print_info : bool = False) -> cupy.ndarray:
    
    # Ref:  https://docs.cupy.dev/en/stable/reference/generated/cupy.linalg.eigh.html
    # note: This function does not allow to select a subset of eig by index.
    #       All the eigenvalues are computed.
    evals, evect = cupy.linalg.eigh( hh )

    mask : cupy.ndarray = (evals == evals[0])
    n_superposition : int = cupy.sum(mask)
    gs : cupy.ndarray = cupy.sum( evect[:,mask].T, axis = 0)/cupy.sqrt(n_superposition)
    
    if (n_superposition == n_eig) and (n_eig != 1):
        print('warning: Number of selected eigenvectors is equal to max allowed eigenvalues. There might be more states with this eigenvalue.')
    if print_info:
        print('[info] gs is superposition of {} states'.format(n_superposition) )
    
    return gs



def timedependent_evolution(htf : callable, T : float, P : int, 
    psi0 : Union[cupy.ndarray, str] = 'superpos', gs_args : dict = {}, 
    force_norm : bool = True, hook_function : Union[callable,None] = None
    ) -> cupy.ndarray:
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
    time = cupy.linspace(0,T,P+1)
    dt = cupy.mean( cupy.diff(time) )
    
    time = cupy.asnumpy( time )
    # fix: this makes the argument of Ht be a numpy value instead of a cupy value
    
    def _single_step_dense(Ht, state):
        evals, evect = cupy.linalg.eigh( Ht )
        U = cupy.exp( -1.j*evals*dt )
        return cupy.dot(evect, U * cupy.dot(evect.transpose().conjugate(), state))
    
    if hook_function is None:
        for t in time[1:]:
            state = _single_step_dense( htf(t), state)
            if force_norm: state = state/(state @ state.conj())
    
    else:
        for t in time[1:]:
            state = _single_step_dense( htf(t), state)
            if force_norm: state = state/(state @ state.conj())
            hook_function(t = t, state = state)
            
    return state


# EXTRAS

def fidelity(psi : cupy.ndarray, chi : cupy.ndarray) -> float:
    cc = psi @ chi.conj()
    fid = cupy.real( cc * cc.conj() )
    return fid.item()