
import numpy as np
import scipy

from typing import Tuple, List, Union

datatype : np.dtype = np.complex128



class operators:
    # Pauli operators
    X : scipy.sparse.spmatrix = scipy.sparse.csr_matrix(np.array([[0, 1], [1, 0]], dtype=datatype))
    Y : scipy.sparse.spmatrix = scipy.sparse.csr_matrix(np.array([[0, -1j], [1j, 0]], dtype=datatype))
    Z : scipy.sparse.spmatrix = scipy.sparse.csr_matrix(np.array([[1, 0], [0, -1]], dtype=datatype))
    
    # trivial operators
    Id : scipy.sparse.spmatrix = scipy.sparse.identity(2, dtype=datatype)
    
    
    
    
class quantum_sparse:
    def __init__(self, size : int, local_dim : int = 2):
        self.size = size
        self.local_dim = local_dim
    
    
    def kron_site(self, operator : scipy.sparse.spmatrix, site : int) -> scipy.sparse.spmatrix:
        if site >= self.size:
            raise Exception('invalid index')
        
        lx : int = site
        rx : int = self.size-site-1
        lo = scipy.sparse.identity(self.local_dim**lx, dtype=datatype) if lx > 0 else 1
        ro = scipy.sparse.identity(self.local_dim**rx, dtype=datatype) if rx > 0 else 1
        
        tmp = scipy.sparse.kron(lo, operator)
        return scipy.sparse.kron(tmp, ro)
    
    
    def kron_list(self, operators : List[scipy.sparse.spmatrix], sites : List[int]) -> scipy.sparse.spmatrix:
        assert len(operators) == len(sites), 'not valid operator sites'
        # sites should be increasingly sorted and without repeated elements
        
        kro = operators[0]
        sites_gap = np.diff( sites ) - 1
        
        for op, gap in zip(operators[1:], sites_gap ):
            if gap > 0:
                kro = scipy.sparse.kron(kro, scipy.sparse.identity(self.local_dim**gap, dtype=datatype) )
            kro = scipy.sparse.kron(kro, op)
        
        # fix L and R offset
        l_offset = sites[0]
        r_offset = self.size - sites[-1] -1
        
        if l_offset != 0:
            kro = scipy.sparse.kron(scipy.sparse.identity(self.local_dim**l_offset, dtype=datatype), kro)
        if r_offset != 0:
            kro = scipy.sparse.kron(kro, scipy.sparse.identity(self.local_dim**r_offset, dtype=datatype))
            
        return kro