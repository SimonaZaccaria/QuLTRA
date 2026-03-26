import numpy as np
import qultra.constants as constants
class J:
    """
    Represents a Josepshon junction component

    Parameters
    ----------
    node_minus : int
        The node to which the negative terminal is connected
    node_plus : int
        The node to which the positive terminal is connected
    J_values: float
        Total linear inductance value [H]. For a junction array, J denotes the inductance of the whole array, not of an individual junction.
    N: int
        Number of junctions. Default is 1.
        (N!=1 implements a JJ array)
    """
    def __init__(self,node_minus,node_plus,J_value,N=1,phi_ext=0):
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.J_value=J_value #[udm]=H
        self.N=N #number of JJ in series
        self.phi_ext=phi_ext
    def admittance(self,z):
        """
        Calculate the admittance of the linear inductor J associated to the junction
        
        Y = 1/s*J

        Parameters
        ----------
        z : complex
            Complex variable (complex frequency)

        Returns
        -------
        Y : complex
            Complex admittance
        """    
        return 1/(z*self.J_value)
    
    def Ej(self):
        """
        Calculate the Josepshon energy of the junction

        Returns
        -------
        Ej : float
            Josephson energy.
        """
        return (constants.phi0/(2*np.pi))**2/self.J_value/constants.h