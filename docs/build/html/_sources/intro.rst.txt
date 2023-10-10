Counterdiabatic Optimized Local Driving annealER
################################################

This package provides tool for Spin system simulation and optimization through a combination of Quantum Optimal Control and Counterdiabatic Driving techniques.
The framework allows to prompt the system through a simplified interface and automatically compute the minimization of the Adiabatic Gauge Potential (AGP) for a given ansatz.



How to install
==============

Clone from GitHub and move inside the project folder

.. code-block:: bash

   git clone https://github.com/CERN-IT-INNOVATION/colder
   cd colder

Install the required packages

.. code-block:: bash

   pip install .




Get started
===========

In this section we will discuss the main interface of `colder`, building examples to gradually introduce the user to its features.
It is strongly suggested to integrate the information of this paragraphs with the provided examples, in order to see the full picture and grasp all the details of this interface.


Prompt systems through hamiltonians
-----------------------------------

The ``physics`` module of `colder` provides the interface to build generic systems of spins. 
The hamiltonian object can be later used for all the next stages of analysis, i.e. building the gauge potential or simulating the quantum annealing schedule.

Let us start from the basics. Import the ``physics`` module.

.. code-block:: python

   import colder.core.physics as cphys

Suppose that you want to build the following hamiltonian,

.. math::

   \mathcal{H}_z = A \cdot \left( \sigma_z^{(0)} + \sigma_z^{(1)} + \sigma_z^{(2)} \right)

which consists of two Pauli Z operators applied on 3 spins, multiplied by a coefficient A. The hamiltonian object would be written in `colder` as

.. code-block:: python

   H_z = cphys.hamiltonian('Z', [0,1,2], coeff='A')


Note that at this stage we are not specifying the total number of spins. This is because the `colder` hamiltonians are just storing the symbols for the later steps.
The user will be asked to provide the total number of spins when needed.

In the same fashion, we can build hamiltonians which involve interaction of more spins. For instance,

.. math::

   \mathcal{H}_{int} = J \cdot \left( \sigma_x^{(0)} \sigma_x^{(1)} + \sigma_x^{(1)} \sigma_x^{(2)} \right)

would be written as

.. code-block:: python

   H_int = cphys.hamiltonian('ZZ', [ (0,1), (1,2) ], coeff='J')

The user can prompt more complex patterns just by changing the disposition of spin indices and the operator strings.
For instance, the hamiltonian

.. math::

   \mathcal{H}_{fancy} = p \cdot \left( \sigma_x^{(0)} \sigma_y^{(1)} + \sigma_y^{(0)} \sigma_x^{(2)} \right)

writes as

.. code-block:: python

   H_fancy = cphys.hamiltonian('XY', [ (0,1), (2,0) ], coeff='p')


Finally, hamiltonian objects can be summed to compose more complex hamiltonians. 
The user can perform this operation either using the ``+`` symbol, or invoking the ``hamiltonian_collection`` class:

.. code-block:: python

   H_sys = H_z + H_int
   # result is the same as
   H_sys = hamiltonian_collection(H_z, H_int)

The user can print on display a dummy expression of the hamiltonian object, using the ``expression()`` method.


Use of lattice module for spin indexing
---------------------------------------

Since the input of spin indices in hamiltonian objects is an essential but tiresome operation, we provide an interface that builds the lists of spins from scratch.
Suppose that you have a 5 spin system, disposed in a linear chain.

.. code-block:: python
   
   >>> import colder.models.lattice as cll
   >>> lattice = cll.chain(nspin = 5)
   >>> lattice.single_site
   [(0,), (1,), (2,), (3,), (4,)]

   >>> lattice.nearest_neighbor
   [(0, 1), (1, 2), (2, 3), (3, 4)]

The module also supports periodic boundary conditions, that can be set using the option ``periodic = True``:

.. code-block:: python
   
   >>> import colder.models.lattice as cll
   >>> lattice = cll.chain(nspin = 5, periodic = True)
   >>> lattice.nearest_neighbor
   [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]


Thus, the input of hamiltonians is simplified:

.. code-block:: python
   
   lattice = cll.chain(nspin = 3, periodic = False)
   H_z = cphys.hamiltonian('Z', lattice.single_site, coeff='A')
   H_int = cphys.hamiltonian('ZZ', lattice.nearest_neighbor, coeff='J')

   H_sys = H_int + H_z


Timedependent coefficient
-------------------------

So far, we have prompted for each term of the hamiltonian a symbolic coefficient through a string.
For instance, in ``H_z = cphys.hamiltonian('Z', lattice.single_site, coeff='A')`` the coefficient is a symbol associated to the letter :math:`A`.
However, in order to execute an annealing schedule, that coefficient has to be a function of time.
In `colder`, the time dependence can be prompted as a regular Python function defined by the user.
The function has to be linked to each term of the hamiltonian using the ``coeff_function`` argument.

As a very minimal example, this is a valid coefficient function initialization:

.. code-block:: python

   tau = 0.1
   def linear_function(t):
      return t/tau
   
   H_test = cphys.hamiltonian('ZZ', lattice.nearest_neighbor,
      coeff = 'J', coeff_function = linear_function
   )

When building more complex systems, there are nevertheless some constraints to follow:
   * The first argument of the function must be the time ``t``. All the following arguments will be treated as extra arguments and have to be passed to the simulation object (see note after the example).
   * All the coefficient functions in a ``hamiltonian_collection`` must have the same I/O.
   * The function has to be vector-safe (i.e. if the input ``t`` is a numpy array, also the return value shall be a numpy array with same dimensions).

.. code-block:: python

   def xf(t, tau):
      return Xfsys*scale_f(t, tau)

   H_X = cphys.hamiltonian('X', lattice.single_site, coeff = 'X', coeff_function=xf)


   def Jf(t, tau):
      return -Jsys*np.ones_like(t)
      #            ^^^ trick to be vector-safe
   def zf(t, tau):
      return Z0sys*np.ones_like(t)
      #            ^^^ trick to be vector-safe

   H_J = cphys.hamiltonian('ZZ', lattice.nearest_neighbor, coeff = 'J', coeff_function=Jf)
   H_Z = cphys.hamiltonian('Z', lattice.single_site, coeff = 'Z', coeff_function=zf)
   
   H = H_J + H_Z + H_X  # the final system hamiltonian

As a final remark, you may notice that all the coefficient functions of ``H`` have the same I/O, and in particular they require an additional argument ``tau``.
When running a `colder` simulation, that extra argument can be prompted using the argument ``system_fargs``, for instance  ``system_fargs = {'tau' : 0.01}``.


Initialize a COLD simulation
----------------------------

After creating the system hamiltonian and the ansatz through the ``physics`` interface, the ``cold`` module allows to initialize a simulation.

.. code-block:: python

   import colder.simulation.cold as ccold

   H = ...
   ansatz = ...

   tau : float = 0.01 # annealing time
   nspin : int = 5

   csym = ccold.cold(
      system = H, ansatz = ansatz, annealing_time = tau, nspin = nspin,
      # the following parameters are to be passed at the system coefficient functions as extra arguments
      system_fargs = {'tau' : tau }
   )

To see the simulation in action, I recommend to take a look at the Ising 1D example notebook.


Complete examples
=================

The user will find Jupyter notebooks in the `examples` folder on the GitHub repository.

* `1D Ising model`_


.. _1D Ising model: https://github.com/CERN-IT-INNOVATION/colder/blob/main/examples/ising-1d.ipynb