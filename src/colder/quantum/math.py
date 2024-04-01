"""This module implements basic functions for quantum information."""

import numpy as np


def fidelity(psi, chi) -> float:
    """
    Returns the fidelity of two quantum states, psi and chi.
    
    .. math:: \mathcal{F} = |\langle\chi|\psi\\rangle|^2

    Parameters
    ----------
    psi: np.ndarray
        quantum state
    chi: np.ndarray
        quantum state

    Returns
    -------
    fid: real
        fidelity of input states
    """
    cc = psi @ chi.conj()
    fid = np.real( cc * cc.conj() )
    return fid

def infidelity(psi, chi) -> float:
    """
    Returns the infidelity of two quantum states, psi and chi.
    
    .. math:: \mathcal{F} = 1 - |\langle\chi|\psi\\rangle|^2

    Parameters
    ----------
    psi: np.ndarray
        quantum state
    chi: np.ndarray
        quantum state

    Returns
    -------
    infid: real
        infidelity of input states
    """
    infid = 1 - fidelity(psi, chi)
    return infid