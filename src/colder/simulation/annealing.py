
import numpy as np

import colder.backend

from typing import Union



class annealing:
    """
    Quantum annealing simulation. The total evolution time is provided by arg `tau`.
    Either one argument must be provided between `P` (total steps) and `dt` (timestep).
    
    Example:
        .. highlight:: python
        .. code-block:: python
        
            import colder.simulation.annealing
            ann = colder.simulation.annealing.annealing(tau = 0.01, P = 200, backend = 'scipy')
            
            def H_t(t):
                # timedependent hamiltonian
                return ...
            
            psi0 = ...  # the initial state
            
            psi_final = ann.run(H_t, psi0=psi0)
    """
    
    
    
    def __init__(self, tau : float, P : Union[int,None] = None, dt : Union[float,None] = None, backend : str = 'scipy', backend_options : dict = {}) -> None:
        """
        Initialize a quantum annealing simulation. The total evolution time is provided by arg `tau`.
        Either one argument must be provided between `P` (total steps) and `dt` (timestep).
        
        Args:
            tau (float): Total annealing time.
            P (Union[int,None], optional): Number of timesteps. Defaults to None. If provided, overrides the value of `dt`.
            dt (Union[float,None], optional): Timestep. Defaults to None. 
            backend (str, optional): Select a backend. Defaults to 'qibo'. Available backends: `qibo`, `scipy`.
            backend_options (dict, optional): Dictionary of kwargs to pass to backend. Defaults to {}.

        Raises:
            Exception: Must provide either the timestep `dt` or the total number of steps `P`. If both values are prompted, `P` is taken in consideration.
            Exception: Not valid backend identifier.
        """
        assert tau > 0, 'total annealing time must be positive'
        if P is None and dt is None:
            raise Exception('must provide at least the number of steps (P) or the timestep (dt).')
        
        self.tau = tau   # total annealing time
        if P is None:
            self.dt = dt     # effective time step
            self.P = tau//dt # number of steps
        else:
            self.P = P       # number of steps
            self.dt = tau/P  # effective time step
        
        self.backend_options = backend_options
        
        if backend == 'qibo':
            self.backend_time_evolution = colder.backend.qibo.routines.timedependent_evolution
        elif backend == 'scipy':
            self.backend_time_evolution = colder.backend.scipy.routines.timedependent_evolution
        #elif backend == 'sparse':
        #    self.backend_time_evolution = colder.backend.sparse.routines.timedependent_evolution
        else:
            raise Exception('not valid backend identifier')
        
    
    
    def run(self, H : callable, psi0 : Union[np.ndarray, str] = 'superpos', hook_function : Union[callable,None] = None) -> np.ndarray:
        """
        Run the annealing simulation with time dependent hamiltonian. H must be a function of time.

        Args:
            H (callable): Function that returns the hamiltonian at time t.
            psi0 (Union[np.ndarray, str], optional): Initial state. Defaults to 'superpos'.
            hook_function (Union[callable,None], optional): Function to execute at each time step. Arguments must be `(t : float, state : np.ndarray)`. Defaults to None.

        Returns:
            np.ndarray: Evolved quantum state.
        """
        return self.backend_time_evolution(H, T = self.tau, P = self.P, psi0 = psi0, *self.backend_options, hook_function = hook_function)