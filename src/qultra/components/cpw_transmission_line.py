import numpy as np
import scipy
import scipy.integrate
import qultra.constants as constants

class CPW:
    """
    Represents a CPW component

    Parameters
    ----------
    node_minus : int
        The node to which the negative terminal is connected
    node_plus : int
        The node to which the positive terminal is connected
    l: float
        Length of the line [m]
    Z0: float
        Charateristic impedence [Ohm]. Default is 50 Ohm
    """
    def __init__(self,node_minus,node_plus,l,Z0=50):
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.l=l #[udm]=m
        self.Z0=Z0 #charaterstic impedence
        self.epsilon_eff=(1+constants.epsilon_r)/2 #for silicon substrate
        self.v=constants.c/np.sqrt(self.epsilon_eff)
        
    
    def admittance_matrix(self,z):
        """
        Construct the complex admittance matrix of a cpw

        Parameters
        ----------
        z : complex
            Complex variable (complex frequency)
        Returns
        -------
        Y_matrix: numpy array    
            complex admittance matrix 
        """
        #Euler exponential
        exp_z=np.exp(self.l*z/self.v)
        exp_mz=np.exp(-self.l*z/self.v)

        #sine and cosine
        sine=(exp_z-exp_mz)/2j
        cosine=(exp_z+exp_mz)/2

        denominator=1j*self.Z0*sine

        #matrix elements
        y_11=cosine/denominator
        y_12=-1/denominator
        y_21=-1/denominator
        y_22=cosine/denominator

        Y_matrix=np.array([[y_11,y_12],[y_21,y_22]])
        return Y_matrix
    
    def current(self,V_0,V_l,z,x):
        """
        Calculate the current in a CPW given the voltage node values.

        Parameters
        ----------
        V_0 : complex
            Voltage at the start node.
        V_l : complex
            Voltage at the end node.
        z : complex
            Complex variable (complex frequency).
        x : float
            Position along the CPW.

        Returns
        -------
        I : complex
            Calculated current at position x.
        """
        exp_z=np.exp(self.l*z/self.v)
        exp_mz=np.exp(-self.l*z/self.v)

        V_minus=(V_l-V_0*exp_mz)/(exp_z-exp_mz)
        V_plus=(V_0*exp_z-V_l)/(exp_z-exp_mz)

        I=(V_plus*np.exp(-x*z/self.v)-V_minus*np.exp(x*z/self.v))/self.Z0
        return I
    
    def inductive_energy(self,V_0,V_l,z):
        """
        Calculate the inductive energy stored in a CPW.

        Parameters
        ----------
        V_0 : complex
            Voltage at the start node.
        V_l : complex
            Voltage at the end node.
        z : complex
            Complex variable (complex frequency).

        Returns
        -------
        E : complex
            Inductive energy stored in the CPW.
        """
        a=0
        b=self.l
        def integrand(x):
            return abs(self.current(V_0,V_l,z,x))**2
        integral_value, _ = scipy.integrate.quad(integrand, a, b) #integrate by using scipy and taking the first value of the tuple
        E=self.Z0*integral_value/(2*self.v)
        return E
