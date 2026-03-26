class C:
    """
    Represent a capacitor component

    Parameters
    ----------
    node_minus : int
        The node to which the negative terminal is connected
    node_plus : int
        The node to which the positive terminal is connected
    C_value : float
        Capacitance value [F]
    """
    def __init__(self,node_minus,node_plus,C_value):
 
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.C_value=C_value #[udm]=F
    
    def admittance(self,z):
        """
        Calculate the admittance of a capacitor C

        Y = z * C

        Parameters
        ----------
        z : complex
            Complex variable (complex frequency)

        Returns
        -------
        Y : complex
            Complex admittance
        """
        return z*self.C_value