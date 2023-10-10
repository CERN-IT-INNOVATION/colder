
import numpy as np
import scipy
import matplotlib.pyplot as plt

from typing import Tuple, List, Union



class interpolator1D:  
    def __init__(self, x_data : np.array, y_data : np.array, bc_type : str = 'not-a-knot'):
        """
        Initialize a 1D spline interpolator from input sample data.

        Args:
            x_data (np.array): Sample data, x.
            y_data (np.array): Sample data, y.
            bc_type (str, optional): Passed to `scipy.interpolate.CubicSpline`. Defaults to 'not-a-knot'.
        """
        # NOTE bc_type = 'clamped' should be optimal for driving terms
        
        self.x_data = x_data
        self.y_data = y_data
        
        assert len(self.x_data) == len(self.y_data), 'input data must be of same len'
        
        self.x_range = tuple([ min(self.x_data), max(self.x_data) ])
        self.n_interp_points = len(self.x_data)
        
        self.spline = scipy.interpolate.CubicSpline(self.x_data, self.y_data, bc_type=bc_type)
        
    @classmethod
    def from_function(cls, function : callable, x_range : tuple, n_interp_points : int = 50, fargs : dict = {}, *args, **kwargs):
        """
        Initialize a 1D spline interpolator from input function and n points sampled in range `x_range`.

        Args:
            function (callable): callable function to interpolate. The arguments must be (`x_data`, fargs).
            x_range (tuple): Tuple containing the range of x data to be sampled.
            n_interp_points (int, optional): How many x data points to sample. Defaults to 50.
            fargs (dict, optional): keyword dictionary passed to `function`. Defaults to {}.

        Returns:
            _type_: interpolator1D
        """
        x_data = np.linspace(x_range[0], x_range[1], n_interp_points)
        y_data = function(x_data, **fargs)
        return cls(x_data, y_data, *args, **kwargs)
    
    
    def get_function(self) -> callable:
        """
        Returns a lambda function calling the spline function. The function first and only argument will be the x point.

        Returns:
            callable: Lambda function wrapping the spline (interpolated) function.
        """
        return lambda x : self.spline(x)
    
    def get_derivative(self, ord : int = 1) -> callable:
        """
        Returns a lambda function calling the spline function for derivative of order `ord`.
        
        Args:
            ord (int, optional): Order of the derivative. Defaults to 1.

        Returns:
            callable: Lambda function wrapping the spline (interpolated) function.
        """
        assert ord > 0, 'must take derivative of order greater than one'
        return lambda x : self.spline(x, ord)
    
    
    
    def plot(self, plot_derivative : bool = True, plot_data : bool = False, labels : List[str] = ["f", "df"], plt_obj = None):
        """
        Plot the interpolated functions.

        Args:
            plot_derivative (bool, optional): If true, the derivatives are plotted. Defaults to True.
            plot_data (bool, optional): If true, plots the original points before interpolation. Defaults to False.
            labels (List[str], optional): Labels for the curve and their derivatives. Defaults to ["f", "df"].
            plt_obj (_type_, optional): Use a matplotlib plt instance instead of the default one. Defaults to None.

        Returns:
            matplotlib.pyplot: `matplotlib.pyplot` figure
        """
        if plt_obj is None:  plt_obj = plt
        
        xs = np.linspace(self.x_range[0], self.x_range[1], self.n_interp_points*20)
        
        plt_obj.plot(xs, self.spline(xs), label=labels[0])
        if plot_derivative:
            plt_obj.plot(xs, self.spline(xs, 1), label=labels[1])
        
        if plot_data:
            plt_obj.plot(self.x_data, self.y_data, 'x', label='spline data')
        
        return plt