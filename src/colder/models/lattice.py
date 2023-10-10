 

# %%
import numpy as np
from typing import Tuple, List, Union


def convert_array_to_target_list(arr : np.array) -> List[tuple]:
    return [ tuple([el]) for el in arr ]

def make_unique_permutations( sites : List[tuple] ) -> List[tuple]:
    pass



class lattice:
    """
    Mother class for any lattice class. Defines the core methods to be inherited by the other classes.
    """
    def __init__(self) -> None:
        """
        Shall be defined specifically for each ineriting class.
        """
        pass
    
    def invert(self, interactions : List[Tuple]) -> List[Tuple]:
        """Inverts list of interaction terms.

        Args:
            interactions (List[Tuple]): List of tuples to invert.

        Returns:
            List[Tuple]: List of inverted interactions.
        """
        return [ tuple(reversed(x)) for x in interactions ]

    def add_invert(self, interactions : List[Tuple]) -> List[Tuple]:
        """Return interaction terms and their inversion.

        Args:
            interactions (List[Tuple]): List of tuples to invert.

        Returns:
            List[Tuple]: List of interactions.
        """
        return interactions + [ tuple(reversed(x)) for x in interactions ]




class chain(lattice):
    """
    1D chain lattice model.
    
    Attributes:
        single_site: single site indices
        nearest_neighbor: nearest neighbor interactions
        next_nearest_neighbor: next nearest neighbor interactions
        
        annni: return both nearest and next nearest neighbor indices
        
        full_connected: pairwise interactions of a fully connected model
        
        three_neighbor: build interactions with 3 adjacent spins
    """
    
    def __init__(self, nspin : int, periodic = False) -> None:
        """
        Initialize a 1D chain with `nspin` spins.

        Args:
            nspin (int): Number of spins.
            periodic (bool, optional): If True, periodic conditions are assumed at boundary. Defaults to False.
        """
        self.nspin = nspin
        
        # all sites
        self.single_site = convert_array_to_target_list( np.arange(nspin) )
        
        # interactions
        self.nearest_neighbor = [ tuple([i,(i+1)%nspin]) for i in range(nspin-1 + periodic) ]
        
        # I'm patching this, can I find a better way?
        if periodic and (nspin < 4):
            self.next_nearest_neighbor = []
        elif periodic and (nspin == 4):
            self.next_nearest_neighbor = [ (0,2), (1,3) ]
        else:
            add = 1 if periodic else -1
            self.next_nearest_neighbor = [ tuple([i-1,(i+1)%nspin]) for i in range(1, nspin+add) ]
            
        #self.next_nearest_neighbor = [ tuple([i,(i+2)%nspin]) for i in range(nspin-2 + periodic) ]
        #self.annni = self.nearest_neighbor + self.next_nearest_neighbor
        
        self.three_neighbor = [ tuple([i,(i+1)%nspin,(i+2)%nspin]) for i in range(nspin-2 + periodic) ]
        self.full_connected = [ tuple([i,j]) for i in range(nspin-1) for j in range(i+1,nspin) ]


class square(lattice):
    """
    2D square lattice model.
    
    Attributes:
        single_site: single site indices
        nearest_neighbor: nearest neighbor interactions
        next_nearest_neighbor: next nearest neighbor interactions
        
        annni: return both nearest and next nearest neighbor indices
        
        full_connected: pairwise interactions of a fully connected model
    """
    
    def __init__(self, shape : tuple) -> None:
        """
        Initialize a 2D square lattice with `nspin` spins.

        Args:
            nspin (int): Number of spins.
            periodic (bool, optional): If True, periodic conditions are assumed at boundary. Defaults to False.
        """
        
        self.nspin = shape[0]*shape[1]
        
        # matrix defining the particle labels
        self.map = np.arange(self.nspin).reshape(shape)
        
        # all sites
        self.single_site = convert_array_to_target_list( np.arange(self.nspin) )
        
        # interactions
        self.horizontal_neighbor = [ (self.map[j,i], self.map[j,i+1]) for i in range(shape[1]-1) for j in range(shape[0]) ]
        self.vertical_neighbor = [ (self.map[j,i], self.map[j+1,i]) for i in range(shape[1]) for j in range(shape[0]-1) ]
        self.nearest_neighbor = self.horizontal_neighbor + self.vertical_neighbor

        self.next_horizontal_neighbor = [ (self.map[j,i], self.map[j,i+2]) for i in range(shape[1]-2) for j in range(shape[0]) ]
        self.next_vertical_neighbor = [ (self.map[j,i], self.map[j+2,i]) for i in range(shape[1]) for j in range(shape[0]-2) ]
        self.next_nearest_neighbor = self.next_horizontal_neighbor + self.next_vertical_neighbor
        
        #self.annni = self.nearest_neighbor + self.next_nearest_neighbor
        self.full_connected = [ tuple([i,j]) for i in range(self.nspin-1) for j in range(i+1,self.nspin) ]
        

