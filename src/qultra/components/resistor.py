class R:
    """
    Represents a resistor component.
 
    Parameters
    ----------
    node_minus : int
        The node to which the negative terminal is connected, the code assumes that at least one of the two nodes is connected to ground.
    node_plus : int
        The node to which the positive terminal is connected, the code assumes that at least one of the two nodes is connected to ground.
    R_value: float
        Resistance value [Ohm]
    """
    def __init__(self,node_minus,node_plus,R_value):
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.R_value=R_value #[udm]=Ohm
    
    def admittance(self,z=None):
        """
        Calculate the admittance of a resistor R
        
        Y = 1/R

        Returns
        -------
        Y : complex
            Complex admittance
        """    
        return 1/self.R_value