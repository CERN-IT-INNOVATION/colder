
import numpy as np


from typing import Tuple, List, Union

datatype : np.dtype = np.complex128



class operators:
    # Pauli operators
    X : np.ndarray = np.array([[0, 1], [1, 0]], dtype=datatype)
    Y : np.ndarray = np.array([[0, -1j], [1j, 0]], dtype=datatype)
    Z : np.ndarray = np.array([[1, 0], [0, -1]], dtype=datatype)
    
    # trivial operators
    Id : np.ndarray = np.identity(2, dtype=datatype)
    
    
    
    
class quantum_sparse:
    def __init__(self, size : int, local_dim : int = 2):
        self.size = size
        self.local_dim = local_dim
    
    
    def kron_site(self, operator : np.ndarray, site : int) -> np.ndarray:
        if site >= self.size:
            raise Exception('invalid index')
        
        lx : int = site
        rx : int = self.size-site-1
        lo = np.identity(self.local_dim**lx, dtype=datatype) if lx > 0 else 1
        ro = np.identity(self.local_dim**rx, dtype=datatype) if rx > 0 else 1
        
        tmp = np.kron(lo, operator)
        return np.kron(tmp, ro)
    
    
    def kron_list(self, operators : List[np.ndarray], sites : List[int]) -> np.ndarray:
        assert len(operators) == len(sites), 'not valid operator sites'
        # sites should be increasingly sorted and without repeated elements
        
        kro = operators[0]
        sites_gap = np.diff( sites ) - 1
        
        for op, gap in zip(operators[1:], sites_gap ):
            if gap > 0:
                kro = np.kron(kro, np.identity(self.local_dim**gap, dtype=datatype) )
            kro = np.kron(kro, op)
        
        # fix L and R offset
        l_offset = sites[0]
        r_offset = self.size - sites[-1] -1
        
        if l_offset != 0:
            kro = np.kron(np.identity(self.local_dim**l_offset, dtype=datatype), kro)
        if r_offset != 0:
            kro = np.kron(kro, np.identity(self.local_dim**r_offset, dtype=datatype))
            
        return kro