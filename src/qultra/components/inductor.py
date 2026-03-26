import numpy as np
import qultra.constants as constants
    
class L:
    """
    Represents an inductor component

    Parameters
    ----------
    node_minus : int
        The node to which the negative terminal is connected
    node_plus : int
        The node to which the positive terminal is connected
    L_value: float
        Inductance value [H]
    """
    def __init__(self,node_minus,node_plus,L_value,phi_ext=0):
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.L_value=L_value #[udm]=H
        self.phi_ext=phi_ext
    
    def admittance(self,z):
        """
        Calculate the admittance of an inductor L

        Y = 1/s*L

        Parameters
        ----------
        z : complex
            Complex variable (complex frequency)

        Returns
        -------
        Y : complex
            Complex admittance
        """    
        return 1/(z*self.L_value)
    
    def El(self):
        """
        Calculate the Inductive energy of the inductor

        Returns
        -------
        El : float
            Inductive energy.
        """
        return (constants.phi0/(2*np.pi))**2/self.L_value/constants.h