
import numpy as np
import scipy
import skopt

import matplotlib.pyplot as plt

import colder.core.physics as cphys
import colder.gauge as cgauge
import colder.backend
import colder.simulation.annealing

from typing import Tuple, List, Union




def get_functions_from_schedule(schedules : dict) -> dict:
    """Extract the function for every interpolation."""
    return { k : v.get_function() for k, v in schedules.items() }


def check_hamiltonian_input(obj) -> Union[cphys.hamiltonian_collection, None]:
    """Check input obj to be hamiltonian_collection or hamiltonian. Otherwise returns None (to be handled)."""
    if isinstance(obj, cphys.hamiltonian_collection):
        return obj
    elif isinstance(obj, cphys.hamiltonian):
        return cphys.hamiltonian_collection(obj)
    else:
        return None





class cold():
    """
    COLD object to handle the annealing simulation and optimization. 
    
    Example:
        .. highlight:: python
        .. code-block:: python
        
            import colder.simulation.cold as ccold
        
            tau : float = 0.01  # annealing time
            nsteps : int = 200  # number of steps for discrete annealing simulation
            Nk : int = 1  # frequencies to control in QOC functions
            
            H = ...
            ansatz = ...
            
            csym = ccold.cold(
                system = H, ansatz = ansatz, annealing_time = tau, nspin = 5,
                system_fargs = {'tau' : tau, 'beta' : np.zeros(Nk)}, 
                # name of the parameter to treat as QOC parameter to optimize
                qoc_parameters = ['beta'],
                backend = 'scipy'
            )
            
            csym.be_ready()
            
            psi0 = csym.get_psi_zero()
            truegs = csym.get_psi_final()

            # define the loss function to optimize
            def loss_infidelity(psi, param):
            # note: the input of this function are the final state psi and the optimized parameter (beta, in this case)
                infidelity = 1 - colder.quantum.math.fidelity(psi, truegs)
                return infidelity
            
            x0 = np.random.uniform(size=Nk) # initial guess of parameters
            optim = csym.run_cold_opti(init_parameters = x0,
                loss_function = loss_infidelity, psi0 = psi0, annealer_nsteps = nsteps,
            )
        
        Remark: the loss function to be optimized must accept as arguments the final state and the optimized parameter.
    """
    
    def __init__(self,
            system : cphys.hamiltonian_collection, ansatz : cphys.hamiltonian_collection,
            annealing_time : float, nspin : int,
            use_cache_agp : Union[bool,cgauge.cache_lcd_numpy,None] = None,
            system_fargs : dict = {}, qoc_parameters : list = None,
            backend : str = 'scipy'
        ) -> None:
        """
        Initialize a COLD annealing simulation.

        Args:
            system (colder.core.physics.hamiltonian_collection): hamiltonian representing the system
            ansatz (colder.core.physics.hamiltonian_collection): hamiltonian representing the ansatz
            annealing_time (float): Total annealing time.
            nspin (int): Number of spins of the system.
            use_cache_agp (Union[bool,cgauge.cache_lcd_numpy,None], optional): If True/None, the adiabatic gauge potential will be cached after its computation. If False, the AGP will not be cached and will be recomputed at every initialization of the simulation. If a `cache_lcd_numpy` is provided, the argument AGP will be explicitly used as cached value. Defaults to None.
            system_fargs (dict, optional): Arguments to pass at system interpolated functions. Defaults to {}.
            qoc_parameters (list, optional): List of string containing the name of the parameter to use as QOC. Defaults to None.
            backend (str, optional): Choose a backend for simulations. Defaults to 'qibo'.
        """
        # processing input args
        self.system : cphys.hamiltonian_collection = check_hamiltonian_input(system)
        if self.system is None:
            raise Exception('system instance is not valid, must be hamiltonian or collection')
        self.system_fargs = system_fargs
        
        self.ansatz : cphys.hamiltonian_collection = check_hamiltonian_input(ansatz)
        if self.ansatz is None:
            raise Exception('ansatz instance is not valid, must be hamiltonian or collection')
        
        self.nspin = nspin
        self.__annealing_time = annealing_time
        
        # manage agp here
        if use_cache_agp is True:
            self.agp = cgauge.adiabaticgp(self.nspin, system, ansatz, enable_cache=False)
        elif use_cache_agp is False:
            self.agp = cgauge.adiabaticgp(self.nspin, system, ansatz, enable_cache=True)
        elif use_cache_agp is None:
            # default class option is used
            self.agp = cgauge.adiabaticgp(self.nspin, system, ansatz)
        elif isinstance(use_cache_agp, cgauge.cache_lcd_numpy):
            self.agp = cgauge.adiabaticgp(self.nspin, system, ansatz, enable_cache=True, inject_cache=use_cache_agp)
        else:
            raise Exception('cache AGP arg is not an adiabatic gauge potential')
        
        
        # manage qoc_parameters
        self.qoc_parameters = qoc_parameters
        if qoc_parameters is not None:
            assert len(qoc_parameters) == 1, 'so far, only one parameter is supported'
            
            # check that system_fargs contains qoc_parameters
            if not all(parname in system_fargs for parname in qoc_parameters):
                raise Exception('missing required control parameter in input qoc list')
        
        # these attributes will be filled with pre-computed results
        self.system_interpolators = None
        self.system_hamiltonians = None
        self.ansatz_interpolators = None
        self.ansatz_hamiltonians = None
        self.__finalize_routine = None
        
        # choose backend
        assert backend in ['scipy', 'qibo'], 'not valid backend'
        if backend == 'qibo':
            self.interface = colder.backend.qibo.interface
            self.routines = colder.backend.qibo.routines
        elif backend == 'scipy':
            self.interface = colder.backend.scipy.interface
            self.routines = colder.backend.scipy.routines
        elif backend == 'sparse': # NOTE: this has been removed in v1.0.2
            self.interface = colder.backend.sparse.interface
            self.routines = colder.backend.sparse.routines
        else:
            raise Exception('not valid backend interface')
        self.backend : str = backend
    
    # function to check runtime requirements
    def __runtime_requirements_int(self) -> bool:
        """
        Check the interpolation & numerical requirements.

        Returns:
            bool: True if the requirements are satisfied.
        """
        requirements = [ self.system_interpolators, self.ansatz_interpolators ]
        return any(elem is None for elem in requirements) 
    
    def __runtime_requirements_ham(self) -> bool:
        """
        Check the hamiltonian requirements.

        Returns:
            bool: True if the requirements are satisfied.
        """
        requirements = [ self.system_hamiltonians, self.ansatz_hamiltonians ]
        return any(elem is None for elem in requirements)
    
    
    def make_schedule_interpolation(self, agp_samples : int = 50) -> None:
        """
        Compute the numerical interpolation of the system driving parameters.

        Args:
            agp_samples (int, optional): Number of samples to use in the interpolation. Defaults to 50.
        """
        
        trange = tuple([0,self.__annealing_time])   # time range of interpolation
        system_interpolators = self.system.make_coefficients_interpolation(trange, self.system_fargs )
        
        finj = self.agp.match_interpolated_injections(system_interpolators) # system functions and their derivatives
        ansatz_interpolators = self.agp.lcd_interpolated_solver(trange=trange, nsamples=agp_samples, finj=finj)
        
        # the interpolators are stored as object attributes
        self.system_interpolators = system_interpolators
        self.ansatz_interpolators = ansatz_interpolators
    
    
    def retrieve_functions_from_schedule(self):
        """Retrieves the function from schedule objects, for both hamiltonian and ansatz."""
        if self.__runtime_requirements_int():
            raise Exception('system_interpolators and/or ansatz_interpolators is None. Run make_schedule_interpolation() first.')
        
        system_coeffs_functions = get_functions_from_schedule(self.system_interpolators)
        ansatz_coeffs_functions = get_functions_from_schedule(self.ansatz_interpolators)
        
        return system_coeffs_functions, ansatz_coeffs_functions
        
    def compute_schedules(self, time : Union[np.ndarray, None], default_time_points : int = 100) -> Union[np.ndarray, np.ndarray]:
        """
        Given an array of time values, computes the schedule function values.

        Args:
            time (Union[np.ndarray,None): Array of time values. If None, the time points will be generated in [0,annealing time].
            default_time_points (int, optional): Number of points to generate, if time array is not provided. Defaults to 100.


        Returns:
            Union[np.ndarray, np.ndarray]: System and ansatz parameters schedule.
        """
        
        if time is None:
            print(f'generating {default_time_points} sample time points')
            time = np.linspace(0, self.__annealing_time, default_time_points)
        
        system_coeffs_functions, ansatz_coeffs_functions = self.retrieve_functions_from_schedule()
        
        sys_sch = np.empty( (len(time), len(system_coeffs_functions)) )
        ans_sch = np.empty( (len(time), len(ansatz_coeffs_functions)) )
        
        for ii, tt in enumerate(time):
            sys_sch[ii, :] = np.array([ ff(tt) for ff in system_coeffs_functions.values() ])
            ans_sch[ii, :] = np.array([ ff(tt) for ff in ansatz_coeffs_functions.values() ])
         
        return sys_sch, ans_sch
    
    
    def plot_schedules(self, npoints = 101, plot_lcd = True, legend = True) -> None:
        """
        Plot the schedule functions.

        Args:
            npoints (int, optional): Number of points to sample. Defaults to 101.
            plot_lcd (bool, optional): If True, plots the LCD terms. Defaults to True.
            legend (bool, optional): If True, legend is enabled. Defaults to True.
        """
        time = np.linspace(0, self.__annealing_time, npoints)
        sys_sch, ans_sch = self.compute_schedules(time)
        
        plt.plot(time, sys_sch, label = self.system_interpolators.keys() )
        if plot_lcd: plt.plot(time, ans_sch, '--', label = self.ansatz_interpolators.keys() )
        
        if legend: plt.legend()
    
    

    def make_hamiltonians(self):
        """Evaluates the static component of hamiltonians (both system and ansatz)."""
        sys_int = self.interface(self.nspin, self.system)
        ans_int = self.interface(self.nspin, self.ansatz)
        
        self.system_hamiltonians = sys_int.COLD_make_hamiltonian()
        self.ansatz_hamiltonians = ans_int.COLD_make_hamiltonian()
        
        self.__finalize_routine = sys_int.COLD_get_finalize_method()
    
    
    def make_timedependent_hamiltonian(self) -> callable:
        """
        Return the timedependent function that computes the hamiltonian for any backend.

        Raises:
            Exception: The requirements to build the hamiltonian are not satisfied. Run `make_hamiltonians()` first.

        Returns:
            callable: Function that return the timedependent hamiltonian at time `t`. Its first and only argument is `t`.
        """
        # collect interpolators and hamiltonians
        system_coeffs_functions, ansatz_coeffs_functions = self.retrieve_functions_from_schedule()
        #   note: this requires the interpolators to be already initialized
        if self.__runtime_requirements_ham():
            raise Exception('system_hamiltonians and/or ansatz_hamiltonians not computed. Run make_hamiltonians() first.')
        
        # note: the .item() in the coefficient is used to get rid of numpy layer
        def tdh(t):
            tsys = sum( system_coeffs_functions[coeff](t).item()*expr for coeff, expr in self.system_hamiltonians )
            tans = sum( ansatz_coeffs_functions[coeff](t).item()*expr for coeff, expr in self.ansatz_hamiltonians )
            
            return self.__finalize_routine(tsys + tans)
        return tdh
    
    
    
    
    def make_timedependent_hamiltonian_nogauge(self) -> callable:
        """
        Return the timedependent function that computes the hamiltonian for any backend. The returned hamiltonian will not have the ansatz AGP terms.

        Raises:
            Exception: The requirements to build the hamiltonian are not satisfied. Run `make_hamiltonians()` first.

        Returns:
            callable: Function that return the timedependent hamiltonian at time `t`. Its first and only argument is `t`.
        """
        # collect interpolators and hamiltonians
        system_coeffs_functions, _ = self.retrieve_functions_from_schedule()
        #   note: this requires the interpolators to be already initialized
        if self.system_hamiltonians is None:
            raise Exception('system_hamiltonians not computed. Run make_hamiltonians() first.')
        
        def tdhg(t):
            terms = sum( system_coeffs_functions[coeff](t).item()*expr for coeff, expr in self.system_hamiltonians ) # FIXME
            return self.__finalize_routine(terms)
        return tdhg
    
    
    def get_psi_zero(self) -> np.ndarray:
        """
        Computes the initial ground state of system hamiltonian via diagonalization.

        Returns:
            np.ndarray: The ground state.
        """
        H0 = self.make_timedependent_hamiltonian_nogauge()
        psi0 = self.routines.get_groundstate_superposition( H0(0) )
        return psi0
    
    def get_psi_final(self) -> np.ndarray:
        """
        Computes the final ground state of system hamiltonian via diagonalization.

        Returns:
            np.ndarray: The ground state (at final annealing time).
        """
        H0 = self.make_timedependent_hamiltonian_nogauge()
        psiF = self.routines.get_groundstate_superposition( H0(self.__annealing_time) )
        return psiF
    
    
    
    
    def get_spectrum(self, n_eig : int, nsteps : int = 100, use_H0 : bool = False) -> Union[np.ndarray, np.ndarray]:
        """
        Returns the spectrum of the system hamiltonian in `nsteps` between time 0 and total annealing time.

        Args:
            n_eig (int): Number of eigenvalues to compute.
            nsteps (int, optional): Number of timesteps. Defaults to 100.
            use_H0 (bool, optional): If True, the hamiltonian is without the AGP. Defaults to False.

        Returns:
            Union[np.ndarray, np.ndarray]: Matrix of eigenvalues and simulation time, respectively.
        """
        if use_H0:
            Ht = self.make_timedependent_hamiltonian_nogauge()
        else:
            Ht = self.make_timedependent_hamiltonian()
            
        time = np.linspace(0,self.__annealing_time, nsteps)
        evalues = []
        for tt in time:
            ev, _ = self.routines.get_lower_eigen(Ht(tt), k = n_eig)
            evalues.append(ev)
            
        return np.matrix(evalues), time
    
    def plot_spectrum(self, n_eig : int = 8, nsteps : int = 100, use_H0 : bool = False) -> None:
        """
        Plots the spectrum of the system hamiltonian in `nsteps` between time 0 and total annealing time.

        Args:
            n_eig (int): Number of eigenvalues to compute.
            nsteps (int, optional): Number of timesteps. Defaults to 100.
            use_H0 (bool, optional): If True, the hamiltonian is without the AGP. Defaults to False.
        """
        ev, time = self.get_spectrum(n_eig = n_eig, nsteps = nsteps, use_H0 = use_H0)
        plt.plot(time, ev)

    
    
    def be_ready(self):
        """Executes all the required pre-computation routines."""
        self.make_schedule_interpolation()
        self.make_hamiltonians()
        
    
    def run(self, nsteps : int = 100, psi0 : Union[np.ndarray, None] = None, dt : Union[float,None] = None, 
            force_be_ready : bool = False, hook_function : Union[callable,None] = None
        ) -> np.ndarray:
        """
        Run an annealing simulation for the current COLD object.

        Args:
            nsteps (int, optional): Number of discretized annealing steps. Defaults to 100.
            psi0 (Union[np.ndarray, None], optional): Initial state. If None, the initial state is compute diagonalizing the hamiltonian without Gauge potential. Defaults to None.
            force_be_ready (bool, optional): If True, be_ready() is executed before the annealing run. Defaults to False.
            hook_function (Union[callable,None], optional): If not None, this function will be called at every step of the annealing schedule. The callable arguments must be `(t : float, state : np.ndarray)`.
        
        Raises:
            Exception: Simulation object is not initialized. Run `be_ready()` method before calling this function to fix this issue.

        Returns:
            np.ndarray: Final annealed state of the system.
        """
        # check be_ready has been executed already (or force it anyways)
        if force_be_ready:
            self.be_ready()
        else:
            if any([self.__runtime_requirements_int(), self.__runtime_requirements_ham()]):
                raise Exception('Schedules and/or hamiltonians not computed. Run be_ready() first.')
        
        if dt is not None:
            # overrides nsteps
            nsteps = None
            
        H_t = self.make_timedependent_hamiltonian()  # timedependent hamiltonian function
        ann = colder.simulation.annealing.annealing(self.__annealing_time, P = nsteps, dt = dt, backend=self.backend) # annealer module
        
        # compute the ground state if not provided as input
        if psi0 is None:   psi0 = self.get_psi_zero()
        
        # execute and return evolved state
        psi = ann.run(H_t, psi0=psi0, hook_function=hook_function)
        return psi
    
    
    
    def make_cold_optimization_function(self, loss_function : callable, psi0 : np.array, annealer_nsteps : int,
            loss_function_args : tuple = ()
        ) -> callable:
        """Creates the function to be optimized, given a specific `loss_function` on the final state.

        Args:
            loss_function (callable): Function to compute the loss. First input is the final state of annealing, second input is the current QOC parameters.
            psi0 (np.array): Initial annealing state.
            annealer_nsteps (int): Number of discrete annealing steps.
            loss_function_args (tuple, optional): Extra arguments to be passed to optimization function. Defaults to ().

        Returns:
            callable: Function to be optimized.
        """
        def f0(params_to_optimize):
            # remake the schedule object with updated parameters
            self.system_fargs[ self.qoc_parameters[0] ] = params_to_optimize
            self.make_schedule_interpolation()
            
            # run the annealing schedule and compute loss function
            psi = self.run(nsteps = annealer_nsteps, psi0=psi0)
            # the loss function args are:  final_annealing_state, parameters, custom_args
            loss = loss_function(psi, params_to_optimize, *loss_function_args)
            return loss
        
        return f0
    
    
    def run_cold_opti(self,
            loss_function : callable, annealer_nsteps : int, psi0 : np.ndarray, 
            init_parameters : np.array, options : dict = {}, method : str = 'Powell',
            loss_function_args : tuple = (), update_param : bool = True
        ) -> scipy.optimize._optimize.OptimizeResult:
        """Run the COLD optimization with backend `scipy.optimize.minimize`.

        Args:
            loss_function (callable): Loss function. Must take two arguments: `psi` (the final state of annealing) and `params` (parameters optimized by the routine).
            annealer_nsteps (int): Number of steps to perform in annealing time evolution.
            psi0 (np.ndarray): Initial state of annealing. It is strongly suggested to compute it manually outside this routine.
            init_parameters (np.array): Initial parameters.
            options (dict, optional): Options to be passed to the optimization routine. Defaults to {}.
            method (str, optional): Method of `scipy.optimize.minimize` to be used. Defaults to 'Powell'.
            loss_function_args (tuple, optional): Extra arguments to be passed to loss_function. Defaults to ().
            update_param (bool, optional): If True, the parameters of the object are update with the result of optimization. Defaults to True.

        Raises:
            Exception: QOC parameters are not specified.
            
        Returns:
            scipy.optimize._optimize.OptimizeResult: Results of optimization, as provided by the backend module `scipy.optimize.minimize`.
        """
        if self.qoc_parameters is None:
            raise Exception('missing qoc parameters, cannot run COLD optimization')
        
        if not update_param:
            # backup the current parameter to be restored later
            bak_system_fargs = self.system_fargs.copy()
        
        f0 = self.make_cold_optimization_function(loss_function=loss_function, psi0=psi0, annealer_nsteps=annealer_nsteps, loss_function_args=loss_function_args)
        
        # deprecated: just Powell optimization, moving to more general minimizer
        #optim = scipy.optimize.fmin_powell(f0, init_parameters, full_output = True, disp = False)
        
        optim = scipy.optimize.minimize(f0, x0=init_parameters, method=method, options=options)
        
        if update_param:
            # update qoc using the optimal parameters
            self.system_fargs[ self.qoc_parameters[0] ] = optim.x
            self.make_schedule_interpolation()
        else:
            # restore old system fargs
            self.system_fargs = bak_system_fargs
            
        return optim



    def run_cold_opti_bayes(self, 
            loss_function : callable, annealer_nsteps : int, psi0 : np.ndarray, 
            dimensions : list, optim_options : dict = {},
            loss_function_args : tuple = (), update_param : bool = True, 
        ) -> scipy.optimize._optimize.OptimizeResult:
        """
        Run the COLD optimization with Bayesian Gaussian backend (`skopt.gp_minimize`).

        Args:
            loss_function (callable): Loss function. Must take two arguments: `psi` (the final state of annealing) and `params` (parameters optimized by the routine).
            annealer_nsteps (int): Number of steps to perform in annealing time evolution.
            psi0 (np.ndarray): Initial state of annealing. It is strongly suggested to compute it manually outside this routine.
            dimensions (list): Bounds on each dimension of the parameters, as required by `skopt.gp_minimize`.
            optim_options (dict, optional): _description_. Defaults to {}. Suggested args to control: `n_calls` and `n_initial_points`.
            loss_function_args (tuple, optional): Extra arguments to be passed to loss_function. Defaults to ().
            update_param (bool, optional): If True, the parameters of the object are update with the result of optimization. Defaults to True.

        Raises:
            Exception: QOC parameters are not specified.

        Returns:
            scipy.optimize._optimize.OptimizeResult: Results of optimization, as provided by the backend module `scipy.optimize.minimize`.
        """
        
        if self.qoc_parameters is None:
            raise Exception('missing qoc parameters, cannot run COLD optimization')
        
        if not update_param:
            # backup the current parameter to be restored later
            bak_system_fargs = self.system_fargs.copy()
        
        f0 = self.make_cold_optimization_function(loss_function=loss_function, psi0=psi0, annealer_nsteps=annealer_nsteps, loss_function_args=loss_function_args)

        # ref: https://scikit-optimize.github.io/stable/modules/generated/skopt.gp_minimize.html
        optim = skopt.gp_minimize(f0,         # the function to minimize
            dimensions = dimensions,    # bounds on each dimension of x
            **optim_options
        )
        
        if update_param:
            # update qoc using the optimal parameters
            self.system_fargs[ self.qoc_parameters[0] ] = optim.x
            self.make_schedule_interpolation()
        else:
            # restore old system fargs
            self.system_fargs = bak_system_fargs
            
        return optim
        