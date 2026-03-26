import numpy as np
import scipy.integrate
import scipy
from scipy.optimize import root
import numbers
from scipy.constants import epsilon_0, mu_0
import qultra.constants as constants



class CPW_coupler:
    """
    Represents a 4-node CPW (coplanar waveguide) coupler.

    This class initializes a CPW coupler component with its node connections,
    gaps, individual CPWs, and physical length. It also computes the 
    capacitance and inductance matrices upon initialization.

    Node ordering (top wiev)
     ::

        1----------2

        3----------4

    The diagrams below illustrate the widths and gaps of the CPW coupler, both in the case without a central ground plane (1) and in the 
    case with a ground plane between the CPW lines (2).

    (1)
     ::

        GND                                           GND
            |       |        |      |        |        |
            |       |        |      |        |        |
            |       |________|      |________|        |  
              gap0    cpw0     gap1    cpw1    gap2


    (2)
     ::

        GND                                                             GND
            |       |        |      |      |        |       |       |
            |       |        |      |      |        |       |       |
            |       |________|      |______|        |_______|       |
              gap0    cpw0     gap1    cpw1   gap2    cpw2    gap3



    Parameters
    ----------
    nodes : list of int
        List of 4 nodes to which the coupler is connected. Must have length 4.
    gap : list of float
        List of gaps between CPWs [um].
    cpw : list
        List of CPW segments' width [um]
    l : float
        Physical length of the coupler [m].

    Attributes
    ----------
    C : ndarray
        Capacitance matrix computed from the CPW configuration.
    L : ndarray
        Inductance matrix computed from the CPW configuration.

    Raises
    ------
    ValueError
        If `nodes` does not contain exactly 4 elements, or if `gap` does not
        have length equal to `len(cpw) + 1`.
    """

    def __init__(self,nodes,gap,cpw,l):
        if len(nodes)!=4:
            raise ValueError ("The component must be connected to 4 nodes")
        if len(gap) != len(cpw) + 1:
            raise ValueError(
                f"`gap` must have length equal to `len(cpw) + 1`. Got len(gap) = {len(gap)}, len(cpw) = {len(cpw)}"
            )
        self.nodes=nodes
        self.gap=gap
        self.cpw=cpw
        self.l=l #[udm]=m
        self.epsilon_eff=(1+constants.epsilon_r)/2 #for silicon substrate
        self.v=constants.c/np.sqrt(self.epsilon_eff)
        self.C, self.L = self.CL_matrices()
        

    
    def branch_point_coordinates(self):
        ''' '''
        gap=self.gap
        cpw=self.cpw
        a=[] #branch points destri
        b=[] #branch points sinistri
        a.append(0) #il punto a_0 ha sempre coordinata 0 che è dove pongo la mia orgine
        b.append(gap[0])
        x=0
        y=gap[0]
        for i in range(len(cpw)):
            x+=(gap[i]+cpw[i])
            y+=(gap[i+1]+cpw[i])
            a.append(x) #coordinate dei punti a che sono alla destra dei conduttori
            b.append(y) #coordinate dei punti b che sono alla sinistra dei conduttori

        a_coordinates=[complex(x) for x in a]
        b_coordinates=[complex(x) for x in b]
        return a_coordinates,b_coordinates

    @staticmethod
    def conformal_mapping(a_coordinates,b_coordinates,c_coordinates):
        def integral_by_part(z0,z1):
            def f1(z):
                u=1
                for c in c_coordinates:
                    u=u*(z-c)
                for a in a_coordinates:
                    if a!=z0 and a!=z1:
                        u=u*(z-a)**(-0.5)
                for b in b_coordinates:
                    if b!=z0 and b!=z1:
                        u=u*(z-b)**(-0.5)
                vf=np.log(z-((z0+z1)/2)+(z-z0)**0.5*(z-z1)**0.5)
                f1=u*vf
                return f1
            def f2(z):
                vf=np.log(z-((z0+z1)/2)+(z-z0)**0.5*(z-z1)**0.5)
                prod1 = np.prod([z - c for c in c_coordinates])
                prod2 = np.prod([(z - a)**(-0.5) for a in a_coordinates if a != z0 and a != z1]) \
                    * np.prod([(z - b)**(-0.5) for b in b_coordinates if b != z0 and b != z1])
                u_prime=0
                for c in c_coordinates:
                    prod3 = np.prod([z - d for d in c_coordinates if d != c])
                    u_prime += prod3 * prod2
                sum_term = 0
                for a in a_coordinates:
                    if a != z0 and a != z1:
                        sum_term += (z - a)**(-1)
                for b in b_coordinates:
                    if b != z0 and b != z1:
                        sum_term += (z - b)**(-1)

                u_prime -= (prod1 * prod2 * sum_term) / 2
                f2=vf*u_prime
                return f2
            integral_part = f1(z1)-f1(z0)
            def z(t):
                return t * (z1 - z0) + z0  # Parametrizzazione del segmento

            real_part, _ = scipy.integrate.quad(lambda t: f2(z(t)).real, 0, 1)
            imag_part, _ = scipy.integrate.quad(lambda t: f2(z(t)).imag, 0, 1)

            numerical_part = real_part + 1j * imag_part
            return integral_part - (z1-z0)*numerical_part   
                        
        ap=[] #nuove coordinate di a dopo il conformal mapping
        bp=[] #nuove coordinate di b dopo il conformal mapping
        ap.append(complex(0)) #il primo resta zero perché è l'origine
        val = integral_by_part(a_coordinates[0], b_coordinates[0])
        bp.append(val)
        for i in range(1,len(b_coordinates)):
            val +=integral_by_part( b_coordinates[i-1], a_coordinates[i])
            ap.append(val)
            val += integral_by_part( a_coordinates[i], b_coordinates[i])
            bp.append(val)
        return ap,bp

    @staticmethod
    def find_c(a_coordinates,b_coordinates, metal_i):
        def c_solve(c_coordinates): #c_coordinates are real number
            def integral_by_part(z0,z1):
                def f1(z):
                    u=1
                    for c in c_coordinates:
                        u=u*(z-complex(c))
                    for a in a_coordinates:
                        if a!=z0 and a!=z1:
                            u=u*(z-a)**(-0.5)
                    for b in b_coordinates:
                        if b!=z0 and b!=z1:
                            u=u*(z-b)**(-0.5)
                    vf=np.log(z-((z0+z1)/2)+(z-z0)**0.5*(z-z1)**0.5)
                    f1=u*vf
                    return f1
                
                def f2(z):
                    vf=np.log(z-((z0+z1)/2)+(z-z0)**0.5*(z-z1)**0.5)
                    prod1 = np.prod([z - complex(c) for c in c_coordinates])
                    prod2 = np.prod([(z - a)**(-0.5) for a in a_coordinates if a != z0 and a != z1]) \
                        * np.prod([(z - b)**(-0.5) for b in b_coordinates if b != z0 and b != z1])
                    u_prime=0
                    for c in c_coordinates:
                        prod3 = np.prod([z - complex(d) for d in c_coordinates if d != c])
                        u_prime += prod3 * prod2
                    sum_term = 0
                    for a in a_coordinates:
                        if a != z0 and a != z1:
                            sum_term += (z - a)**(-1)
                    for b in b_coordinates:
                        if b != z0 and b != z1:
                            sum_term += (z - b)**(-1)

                    u_prime -= (prod1 * prod2 * sum_term) / 2
                    f2=vf*u_prime
                    return f2
                
                integral_part = f1(z1)-f1(z0)
                def z(t):
                    return t * (z1 - z0) + z0  # Parametrizzazione del segmento

                real_part, _ = scipy.integrate.quad(lambda t: f2(z(t)).real, 0, 1)
                imag_part, _ = scipy.integrate.quad(lambda t: f2(z(t)).imag, 0, 1)

                numerical_part = real_part + 1j * imag_part
                return integral_part - (z1-z0)*numerical_part
            
            constraints=[]
            for j in range(len(a_coordinates)):
                if j!=metal_i and j!=(metal_i+1):
                    f=integral_by_part(a_coordinates[j],b_coordinates[j])
                    constraints.append(f.imag)
            return constraints
        
        c_coordinates_guest=[]
        for j in range(len(a_coordinates)):
            if j!=metal_i and j!=(metal_i+1):
                c=(b_coordinates[j].real+a_coordinates[j].real)/2
                c_coordinates_guest.append(c)
        sol=root(c_solve,c_coordinates_guest)
        if not sol.success:
            raise RuntimeError("Root finding failed: " + sol.message)
        c_coordinates=[complex(x) for x in sol.x]
        return c_coordinates           

    
    def CL_matrices(self):
        ''' '''
        gap=self.gap
        cpw=self.cpw
        a,b=self.branch_point_coordinates()
        C=np.zeros((len(cpw),len(cpw)))
        for j in range(len(cpw)):
            c=self.find_c(a,b,j)
            ap,bp=self.conformal_mapping(a,b,c)
            for i in range(len(cpw)):
                C[i,j]=(constants.epsilon_r+1)*constants.epsilon_0*(bp[i].real-ap[i+1].real)/bp[j].imag
        if C.shape[0] == 3:
            C = np.delete(C, 1, axis=0)
            C = np.delete(C, 1, axis=1)

        L=np.linalg.inv(C)/self.v**2

        return C,L

    def Y(self,z):
        """
        Construct the complex admittance matrix of a cpw coupler

        Parameters
        ----------
        z : complex
            Complex variable (complex frequency)
        Returns
        -------
        Y_matrix: numpy array    
            complex admittance matrix 
        """

        l = self.l
        C= self.C

        Z_matrix_inv = self.v* C
        exp_z = np.exp(l * z / self.v)
        exp_mz = np.exp(-l * z / self.v)
        imag_part = exp_z - exp_mz
        dim = Z_matrix_inv.shape[0]

        Y_matrix = np.zeros((2 * dim, 2 * dim), dtype=np.complex128)

        for i in range(2 * dim):
            k = i // 2   # indice della cella
            even = (i % 2 == 0)

            if even:
                V_plus_k = exp_z / imag_part
                V_minus_k = -exp_mz / imag_part
            else:
                V_plus_k = -1 / imag_part
                V_minus_k = 1 / imag_part

            # V_plus and V_minus are zero everywhere except at index k
            I_plus = Z_matrix_inv[:, k] * V_plus_k
            I_minus = -Z_matrix_inv[:, k] * V_minus_k

            for j in range(dim):
                Y_matrix[2 * j, i] = I_plus[j] + I_minus[j]
                Y_matrix[2 * j + 1, i] = -(I_plus[j] * exp_mz + I_minus[j] * exp_z)

        return Y_matrix
    
    
    def inductive_energy(self,V,z):
        """
        Calculate the inductive energy stored in a CPW.

        Parameters
        ----------

        V : complex array
            Voltage at nodes.
        z : complex
            Complex variable (complex frequency).

        Returns
        -------
        E : complex
            Inductive energy stored in the CPW.

        """
        a=0
        b=self.l
        L=self.L   

        def integrand(x):
            I=self.current(V,z,x)
            dW=(I.conj().T @ L @ I)/2
            return np.real(dW.item())
        E, _ = scipy.integrate.quad(integrand, a, b) #integrate by using scipy and taking the first value of the tuple
        return E
    
    def current(self,V,z,x):
            
        """
        Calculate the current in a CPW coupler given the voltage node values.

        Parameters
        ----------
        V : complex array
            Voltage at  nodes.
        z : complex
            Complex variable (complex frequency).
        x : float
            Position along the CPW coupler.

        Returns
        -------
        I : complex
            Calculated current at position x.
        """
        C = self.C

        Z_matrix_inv = self.v* C
        dim = Z_matrix_inv.shape[0]
        exp_z=np.exp(self.l*z/self.v)
        exp_mz=np.exp(-self.l*z/self.v)
        imag_part=exp_z-exp_mz

        V_plus=np.zeros((dim,1),dtype=np.complex128)
        V_minus=np.zeros((dim,1),dtype=np.complex128)

        for j in range(dim):
            V_plus[j,0]=(V[2*j]*exp_z-V[2*j+1])/imag_part
            V_minus[j,0]=(-V[2*j]*exp_mz+V[2*j+1])/imag_part

        I_plus=Z_matrix_inv @ V_plus
        I_minus=-Z_matrix_inv @ V_minus

        I = I_plus * np.exp(-x * z / self.v) + I_minus * np.exp(x * z / self.v)

        return I    