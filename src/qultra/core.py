import numpy as np
import scipy
import scipy.integrate
import numbers
from tabulate import tabulate

try:
    from .constants import *   # import relativo (da pacchetto)
except ImportError:
    from constants import *   # import assoluto (da sorgente)


try:
    from .find_zeros import*
except ImportError:
    # When running from source without pip installation
    from find_zeros import *
'''
def zero_algo(admittance, starting_point, final_point):
     
        Finds a zero (root) of a complex admittance function in the jω domain.

        This function takes a complex-valued admittance function defined in the jω domain 
        and attempts to find a frequency ω such that admittance(jω) = 0 within the interval 
        [starting_point, final_point]. It is typically used to compute the resonant or 
        eigenfrequencies of a circuit by passing the characteristic polynomial as input.

        Args:
            admittance: A callable function of a complex variable z = jω.
            starting_point:The lower bound of the frequency range (in GHz).
            final_point: The upper bound of the frequency range (in GHz).

        Returns:
            The frequency f [GHz] where admittance(jω) ≈ 0, if found.

        Raises:
            ValueError: If root-finding fails or no zero is detected in the interval.
     

     minimum_points=[]
     step_for_interval=0.1
    
     for start_interval in np.arange(starting_point,final_point,step_for_interval):
    
        step_for_nodes=0.001
        minimum_value=float('inf')
        end_interval=start_interval+step_for_interval+step_for_nodes
        N=int((end_interval-start_interval)/step_for_nodes)+1
    
        for f in np.linspace(start_interval,end_interval,N):
            function_value=abs(admittance(1j*2*np.pi*1e9*f))
            if function_value< minimum_value:
                minimum_value=function_value
                f_min=f
        if f_min != start_interval and f_min !=end_interval:
            
            new_start=f_min-step_for_nodes
            new_end=f_min+step_for_nodes
            new_step_for_nodes=step_for_nodes/10
            k=1
            N=int((new_end-new_start)/new_step_for_nodes)+1
            minimum_value=float('inf')
            minimum_value_array=[]
            while k<26:
                for f in np.linspace(new_start,new_end,N):
                    function_value=abs(admittance(1j*2*np.pi*1e9*f))
                    if function_value< minimum_value:
                        minimum_value=function_value
                        f_min=f
                new_start=f_min-new_step_for_nodes
                new_end=f_min+new_step_for_nodes
                new_step_for_nodes=new_step_for_nodes/10
                N=int((new_end-new_start)/new_step_for_nodes)+1
                minimum_value_array.append(minimum_value)
                k=k+1
            if minimum_value_array[24]*1e5<minimum_value_array[0]:
                minimum_points.append(f_min)
     if not minimum_points:
         raise ValueError("No zeros found in the specified interval.")
    
     return minimum_points

def zero_algo_complete(admittance,f_starting_point, f_end_point):
    """
    Finds complex roots (zeros) of a complex admittance function within a given frequency range.

    This function searches for local minima of |admittance(z)| where z = k + jω,
    sweeping through a specified real frequency interval. For each small interval,
    it refines the search in both real and imaginary components to identify precise
    complex frequencies where the magnitude of the admittance approaches zero.

    Args:
        admittance: A callable function of a complex variable z, representing the admittance of the system.
        f_starting_point: The lower bound of the frequency range (in GHz).
        f_end_point: The upper bound of the frequency range (in GHz).

    Returns:
        A list of complex frequencies [f,k] where the admittance function is (approximately) zero.
        f is expressed in GHz, k is expressed in MHz.
    
    Raises:
        ValueError: If no zero is found within the specified frequency range.
    """

    minimum_points=[]
    f_step=0.01
    k_step=0.01

    k_starting_point=-10
    k_end_point=1 #MHz

    N_f=int((f_end_point-f_starting_point)/(100*f_step))+1
    N_k=int((k_end_point-k_starting_point)/(100*k_step))+1
    

    for f_start_interval in np.linspace(f_starting_point,f_end_point,N_f):
        
         f_end_interval=f_start_interval+101*f_step
         for k_start_interval in np.linspace(k_starting_point,k_end_point,N_k):
             k_end_interval=k_start_interval+101*k_step
             minimum_value=float('inf')
             N1_f=int((f_end_interval-f_start_interval)/f_step)+1
             N1_k=int((k_end_interval-k_start_interval)/k_step)+1
             for f in np.linspace(f_start_interval,f_end_interval,N1_f):
    
                 for k in np.linspace(k_start_interval,k_end_interval,N1_k):
                     z=2*np.pi*1e6*k+1j*2*np.pi*1e9*f #f in GHz, k in MHz
                     function_value=abs(admittance(z))
                
                     if function_value< minimum_value:
                        minimum_value=function_value
                        f_min=f
                        k_min=k
            
             if f_min != f_start_interval and f_min !=f_end_interval and k_min !=k_start_interval and k_min!=k_end_interval:
                 f_new_start=f_min-f_step
                 f_new_end=f_min+f_step
                 k_new_start=k_min-k_step
                 k_new_end=k_min+k_step
                 f_new_step=f_step/10
                 k_new_step=k_step/10
                 j=1
                 minimum_value=float('inf')
                 minimum_value_array=[]
                 N2_f=int((f_new_end-f_new_start)/f_new_step)+1
                 N2_k=int((k_new_end-k_new_start)/k_new_step)+1
                 while j<26:
                     for f in np.linspace(f_new_start,f_new_end,N2_f):
                         for k in np.linspace(k_new_start,k_new_end,N2_k):
                             z=2*np.pi*1e6*k+1j*2*np.pi*1e9*f #f in GHz, k in MHz
                             function_value=abs(admittance(z))
                             if function_value< minimum_value:
                                 minimum_value=function_value
                                 f_min=f
                                 k_min=k
                     f_new_start=f_min-f_new_step
                     f_new_end=f_min+f_new_step
                     k_new_start=k_min-k_new_step
                     k_new_end=k_min+k_new_step
                     f_new_step=f_new_step/10
                     k_new_step=k_new_step/10
                     minimum_value_array.append(minimum_value)
                     j=j+1
                 
                 if minimum_value_array[24]*1e3<minimum_value_array[0]:
                     minimum_points.append([f_min,k_min])
    
    if not minimum_points:
        raise ValueError("No zeros found in the specified interval.")
    return minimum_points

def zero_algo(admittance, starting_point, final_point):
    
    Finds a zero (root) of a complex admittance function in the jω domain.

    This function takes a complex-valued admittance function defined in the jω domain 
    and attempts to find a frequency ω such that admittance(jω) = 0 within the interval 
    [starting_point, final_point]. It is typically used to compute the resonant or 
    eigenfrequencies of a circuit by passing the characteristic polynomial as input.

    Args:
        admittance: A callable function of a complex variable z = jω.
        starting_point:The lower bound of the frequency range (in GHz).
        final_point: The upper bound of the frequency range (in GHz).

    Returns:
        The frequency f [GHz] where admittance(jω) ≈ 0, if found.

    Raises:
        ValueError: If root-finding fails or no zero is detected in the interval.
    
    def f(z):
        return admittance(2*np.pi*1e9*z)
    R=cxroots.Rectangle([-0.5,0.5],[starting_point, final_point])
    sol=R.roots(f)
    minimum_points=[zero.imag for zero in sol.roots]
    return minimum_points

def zero_algo_complete(admittance,starting_point, final_point):
    def f(z):
        return admittance(2*np.pi*1e6*z.real+2j*np.pi*1e9*z.imag)
    R=cxroots.Rectangle([-10/1e3,0],[starting_point, final_point])
    sol=R.roots(f,tol=1e-5)
    minimum_points=[[zero.imag,zero.real] for zero in sol.roots]
    return minimum_points
'''
class C:
    def __init__(self,node_minus,node_plus,C_value):
        """
        Represents a capacitor component

        Parameters:
        node_minus, node_plus: The nodes to which the component is connected
        C_value: Capacitance value [F]
        """
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.C_value=C_value #[udm]=F
    
    def admittance(self,z):
        """
        Calculate the admittance of a capacitor C
        Y=s*C
        Args:
            s:complex variable
        Returns:
            complex admittance Y
        """    
        return z*self.C_value
    
class L:
    """
    Represents an inductor component

    Parameters:
    node_minus, node_plus: The nodes to which the component is connected
    L_value: Inductance value [H]
    """
    def __init__(self,node_minus,node_plus,L_value):
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.L_value=L_value #[udm]=H
    
    def admittance(self,z):
        """
        Calculate the admittance of an inductor L
        Y=1/s*L
        Args:
            s:complex variable
        Returns:
            complex admittance Y
        """    
        return 1/(z*self.L_value)
    
class R:
    """
    Represents a resistor component

    Parameters:
    node_minus, node_plus: The nodes to which the component is connected
    R_value: Resistance value [Ohm]
    """
    def __init__(self,node_minus,node_plus,R_value):
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.R_value=R_value #[udm]=Ohm
    
    def admittance(self,z=None):
        """
        Calculate the admittance of a resistor R
        Y=1/R
        Returns:
            complex admittance Y
        """    
        return 1/self.R_value
class J:
    """
    Represents a Josepshon junction component

    Parameters:
    node_minus, node_plus: The nodes to which the component is connected
    J_values: LInear inductance value [H]
    N: Number of junctions (N!=1 implements a JJ array)
    """
    def __init__(self,node_minus,node_plus,J_value,N=1):
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.J_value=J_value #[udm]=H
        self.N=N #number of JJ in series
    
    def admittance(self,z):
        """
        Calculate the admittance of an inductor L
        Y=1/s*L
        Args:
            z:complex variable
        Returns:
            complex admittance Y
        """    
        return 1/(z*self.J_value)
    
    def Ej(self):
        """
        Calculate the Josepshon energy of the junction
        """
        return (phi0/(2*np.pi))**2/self.J_value/h

class CPW:
    """
    Represents a CPW component

    Parameters:
    node_minus, node_plus: The nodes to which the component is connected
    l: Length of the line
    Z0: Charateristic impedence
    """
    def __init__(self,node_minus,node_plus,l,Z0=50):
        self.node_minus=node_minus
        self.node_plus=node_plus
        self.l=l #[udm]=m
        self.Z0=Z0 #charaterstic impedence
        
    
    def admittance_matrix(self,z):
        """
        Construct the complex admittance matrix of a cpw

        Args:
            z: complex variable
        Returns:
            complex admittance matrix 
        """
        #Euler exponential
        exp_z=np.exp(self.l*z/v)
        exp_mz=np.exp(-self.l*z/v)

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
        Calculate current in a cpw given the voltage node values
        """
        exp_z=np.exp(self.l*z/v)
        exp_mz=np.exp(-self.l*z/v)

        V_minus=(V_l-V_0*exp_mz)/(exp_z-exp_mz)
        V_plus=(V_0*exp_z-V_l)/(exp_z-exp_mz)

        I=(V_plus*np.exp(-x*z/v)-V_minus*np.exp(x*z/v))/self.Z0
        return I
    
    def inductive_energy(self,V_0,V_l,z):
        """
        Calculate the inductive energy stored in a cpw
        """
        a=0
        b=self.l
        def integrand(x):
            return abs(self.current(V_0,V_l,z,x))**2
        integral_value, _ = scipy.integrate.quad(integrand, a, b) #integrate by using scipy and taking the first value of the tuple
        E=self.Z0*integral_value/(2*v)
        return E


class CPW_coupler:
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
        self.C, self.L = self.CL_matrices()

    def branch_point_coordinates(self):
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
                v=np.log(z-((z0+z1)/2)+(z-z0)**0.5*(z-z1)**0.5)
                f1=u*v
                return f1
            def f2(z):
                v=np.log(z-((z0+z1)/2)+(z-z0)**0.5*(z-z1)**0.5)
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
                f2=v*u_prime
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
                    v=np.log(z-((z0+z1)/2)+(z-z0)**0.5*(z-z1)**0.5)
                    f1=u*v
                    return f1
                
                def f2(z):
                    v=np.log(z-((z0+z1)/2)+(z-z0)**0.5*(z-z1)**0.5)
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
                    f2=v*u_prime
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
        gap=self.gap
        cpw=self.cpw
        a,b=self.branch_point_coordinates()
        C=np.zeros((len(cpw),len(cpw)))
        for j in range(len(cpw)):
            c=self.find_c(a,b,j)
            ap,bp=self.conformal_mapping(a,b,c)
            for i in range(len(cpw)):
                C[i,j]=(epsilon_r+1)*epsilon_0*(bp[i].real-ap[i+1].real)/bp[j].imag
        if C.shape[0] == 3:
            C = np.delete(C, 1, axis=0)
            C = np.delete(C, 1, axis=1)

        L=np.linalg.inv(C)/v**2

        return C,L

    def Y(self,z):
        l = self.l
        C= self.C

        Z_matrix_inv = v* C
        exp_z = np.exp(l * z / v)
        exp_mz = np.exp(-l * z / v)
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
            Calculate the inductive energy stored in a cpw
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
            Calculate current in a cpw given the voltage node values
            """
            C = self.C

            Z_matrix_inv = v* C
            dim = Z_matrix_inv.shape[0]
            exp_z=np.exp(self.l*z/v)
            exp_mz=np.exp(-self.l*z/v)
            imag_part=exp_z-exp_mz

            V_plus=np.zeros((dim,1),dtype=np.complex128)
            V_minus=np.zeros((dim,1),dtype=np.complex128)

            for j in range(dim):
                V_plus[j,0]=(V[2*j]*exp_z-V[2*j+1])/imag_part
                V_minus[j,0]=(-V[2*j]*exp_mz+V[2*j+1])/imag_part

            I_plus=Z_matrix_inv @ V_plus
            I_minus=-Z_matrix_inv @ V_minus

            I = I_plus * np.exp(-x * z / v) + I_minus * np.exp(x * z / v)

            return I      

class QCircuit:
    """
    Represents a quantum circuit

    Parameters:
    netlist: the list of component instances that define the circuit
    f_starting_point, f_end_point: the start and end frequencies of the interval over which the circuit is analyzed
    """

    def __init__(self, netlist, f_starting_point,f_end_point):
        self.netlist=netlist
        self.f_starting_point = f_starting_point
        self.f_end_point = f_end_point

        if len(self.netlist) == 0:
            raise ValueError("There are no components in the circuit")
        if self.shorts():
            raise ValueError("Your circuit appears to be open or shorted making the analysis impossible")
        if not self.is_connected():
            raise ValueError("Your circuit appears to be not connected making the analysis impossible")

        self.modes=self.eigenvalues(self.f_starting_point, self.f_end_point)
        

    def shorts(self):
        """
        Check if the circuit is shorted
        """
        components=self.netlist
        for comp in components:
            if hasattr(comp,"node_minus"):
                if comp.node_minus == comp.node_plus:
                    return True
            if hasattr(comp, "nodes"):
                if (comp.nodes[0]==comp.nodes[1]) and (comp.nodes[2]== comp.nodes[3]) and (comp.nodes[0] == comp.nodes[2]):
                    return True
            
        return False
    '''
    def is_connected(self):
        """
        Check if the circuit is connected
        """
    # Collect all nodes present in the components
        components=self.netlist
        all_nodes = []
        for comp in components:
            if comp.node_minus not in all_nodes:
                all_nodes.append(comp.node_minus)
            if comp.node_plus not in all_nodes:
                all_nodes.append(comp.node_plus)

        # Check if ground node (0) is present
        if 0 not in all_nodes:
            return False  # No ground node means circuit is not connected

        # Build adjacency list: for each node, list its connected nodes
        connections = {}
        for comp in components:
            n1 = comp.node_minus
            n2 = comp.node_plus
            if n1 not in connections:
                connections[n1] = []
            if n2 not in connections:
                connections[n2] = []
            connections[n1].append(n2)
            connections[n2].append(n1)

        # Initialize lists for nodes visited and nodes to visit
        visited = []
        to_visit = [0]  # Start traversal from ground node (0)

        # Perform simple depth-first search (DFS) to find reachable nodes
        while to_visit:
            current = to_visit.pop()  # Take a node to visit
            if current not in visited:
                visited.append(current)  # Mark node as visited
                # Add all adjacent nodes not yet visited to to_visit list
                for neighbor in connections.get(current, []):
                    if neighbor not in visited:
                        to_visit.append(neighbor)

        # Check if all nodes were visited (reachable from ground)
        for node in all_nodes:
            if node not in visited:
                return False  # Node not reachable, circuit not fully connected

        return True  # All nodes reachable, circuit is connected
    '''

    def is_connected(self):
        """
        Check if the circuit is connected
        """
        components = self.netlist

        # Collect all nodes present in the components
        all_nodes = []
        for comp in components:
            if hasattr(comp, "node_minus") and hasattr(comp, "node_plus"):
                if comp.node_minus not in all_nodes:
                    all_nodes.append(comp.node_minus)
                if comp.node_plus not in all_nodes:
                    all_nodes.append(comp.node_plus)
            elif hasattr(comp, "nodes"):
                for node in comp.nodes:
                    if node not in all_nodes:
                        all_nodes.append(node)

        # Check if ground node (0) is present
        if 0 not in all_nodes:
            return False  # No ground node means circuit is not connected

        # Build adjacency list: for each node, list its connected nodes
        connections = {}
        for comp in components:
            if hasattr(comp, "node_minus") and hasattr(comp, "node_plus"):
                n1 = comp.node_minus
                n2 = comp.node_plus
                if n1 not in connections:
                    connections[n1] = []
                if n2 not in connections:
                    connections[n2] = []
                connections[n1].append(n2)
                connections[n2].append(n1)
            elif hasattr(comp, "nodes"):
                nodes = comp.nodes
                # Connect nodes in pairs (assume CPW coupler connects 0-1 and 2-3)
                if len(nodes) == 4:
                    pairs = [(0, 1), (2, 3)]
                    for i, j in pairs:
                        n1 = nodes[i]
                        n2 = nodes[j]
                        if n1 not in connections:
                            connections[n1] = []
                        if n2 not in connections:
                            connections[n2] = []
                        connections[n1].append(n2)
                        connections[n2].append(n1)

        # Initialize lists for nodes visited and nodes to visit
        visited = []
        to_visit = [0]  # Start traversal from ground node (0)

        # Perform simple depth-first search (DFS) to find reachable nodes
        while to_visit:
            current = to_visit.pop()
            if current not in visited:
                visited.append(current)
                for neighbor in connections.get(current, []):
                    if neighbor not in visited:
                        to_visit.append(neighbor)

        # Check if all nodes were visited (reachable from ground)
        for node in all_nodes:
            if node not in visited:
                return False  # Node not reachable, circuit not fully connected

        return True  # All nodes reachable, circuit is connected

    def build_total_Y_matrix(self, z):
        """
        Builds the total admittance matrix Y for the circuit,
        assuming node 0 is ground and should be excluded from the final matrix.
        
        Args:
            components: list of component instances (C, L, R, J, CPW)
            z: complex variable (e.g., s = jω or k + jω)
        
        Returns:
            Reduced admittance matrix (excluding ground node 0)
        """
        components=self.netlist
        # Determine the highest node number
        #max_node = 0
        #for comp in components:
        #    max_node = max(max_node, comp.node_minus, comp.node_plus)
        max_node = 0
        for comp in components:
            if hasattr(comp, "node_minus") and hasattr(comp, "node_plus"):
                max_node = max(max_node, comp.node_minus, comp.node_plus)
            elif hasattr(comp, "nodes"):
                max_node = max(max_node, *comp.nodes)


        N_total = max_node+1  # total number of nodes including ground
        Y_total = np.zeros((N_total, N_total), dtype=complex)

        for comp in components:
            if hasattr(comp, "node_minus") and hasattr(comp, "node_plus"):
                n1 = comp.node_minus
                n2 = comp.node_plus
            """"
            # Get admittance or local Y matrix
            if isinstance(comp, C):
                Y = comp.C_admittance(z)
            elif isinstance(comp, L):
                Y = comp.L_admittance(z)
            elif isinstance(comp, R):
                Y = comp.R_admittance()
            elif isinstance(comp, J):
                Y = comp.J_admittance(z)
            elif isinstance(comp, CPW):
                Y_local = comp.CPW_admittance_matrix(z)
                Y_total[n1, n1] += Y_local[0, 0]
                Y_total[n1, n2] += Y_local[0, 1]
                Y_total[n2, n1] += Y_local[1, 0]
                Y_total[n2, n2] += Y_local[1, 1]
                continue
            else:
                raise TypeError(f"Unsupported component type: {type(comp)}")

            # Fill the global Y matrix (only for scalar admittances)
            if n1 != n2:
                Y_total[n1, n1] += Y
                Y_total[n2, n2] += Y
                Y_total[n1, n2] -= Y
                Y_total[n2, n1] -= Y
            else:
                Y_total[n1, n1] += Y

            """
            if hasattr(comp, "admittance_matrix"):
                Y_matrix = comp.admittance_matrix(z)
                Y_total[n1, n1] += Y_matrix[0, 0]
                Y_total[n1, n2] += Y_matrix[0, 1]
                Y_total[n2, n1] += Y_matrix[1, 0]
                Y_total[n2, n2] += Y_matrix[1, 1]
            elif hasattr(comp, "admittance"):
                Y = comp.admittance(z)
                
                Y_total[n1, n1] += Y
                Y_total[n2, n2] += Y
                Y_total[n1, n2] -= Y
                Y_total[n2, n1] -= Y
            
            
            elif hasattr(comp,"Y"):
                Y_matrix=comp.Y(z)
                n=len(comp.nodes)
                for i in range(n):
                    for j in range(n):
                        Y_total[comp.nodes[i],comp.nodes[j]]+=Y_matrix[i,j]

                
            else:
                raise TypeError(f"Component {comp} has no admittance method or admittance matrix object")
        # Remove row and column corresponding to ground node (node 0)
        Y_reduced = Y_total[1:, 1:] #to eliminate ground row and coulomn
       
        return Y_reduced

    def characteristic_polynomial_reduced(self,z):
        nodes_to_delete=[]
        for comp in self.netlist:
            if isinstance(comp,R):
                nodes_to_delete.append(max(comp.node_minus,comp.node_plus)-1)

        Y_matrix=self.build_total_Y_matrix(z)
        Y_matrix=np.delete(Y_matrix, nodes_to_delete, axis=0)
        Y_matrix=np.delete(Y_matrix, nodes_to_delete, axis=1)
        det_Y=np.linalg.det(Y_matrix)

        K = scipy.linalg.null_space(Y_matrix,rcond=1e-10)
        dim_kernel = K.shape[1]
        
        return det_Y, Y_matrix.shape[0], dim_kernel
    '''
    def check_singularities(self, nodes_to_delete):
        ranks = []
        f_list = np.arange(1, 2.1, 0.1)
        n= None
        for f in f_list:
            Y_matrix = self.build_total_Y_matrix(2j*np.pi*1e9*f)
            Y_matrix=np.delete(Y_matrix, nodes_to_delete, axis=0)
            Y_matrix=np.delete(Y_matrix, nodes_to_delete, axis=1)

            if not np.isfinite(Y_matrix).all():
                return True
            if n is None:
                n = Y_matrix.shape[0]  # Salvo la dimensione la prima volta

            rank = np.linalg.matrix_rank(Y_matrix)
            print(rank)
            print(np.linalg.det(Y_matrix))
            if rank==n:
                return False
        return True
    '''
    def characteristic_polynomial(self,z):
        Y_matrix=self.build_total_Y_matrix(z)
        det_Y=np.linalg.det(Y_matrix)
        K = scipy.linalg.null_space(Y_matrix,rcond=1e-10)
        dim_kernel = K.shape[1]
        return det_Y, Y_matrix.shape[0],dim_kernel
    
    def there_is_R(self):
        """
        Check if there are resistive components
        """
        for comp in self.netlist:
            if isinstance(comp, R):
                return True
            
        return False
    
    def eigenvalues(self,f_starting_point, f_end_point):
        """
        Find the eigenfrequencies (modes) of the circuit
        """

        if f_starting_point >= f_end_point:
            raise ValueError("f_starting_point must be < f_end_point")
        if f_starting_point<0 or f_end_point<0:
             raise ValueError("frequencies must be positive")
    
        if self.there_is_R():
            guesses=zero_algo(self.characteristic_polynomial_reduced,f_starting_point, f_end_point)
            modes=zero_algo_complete(self.characteristic_polynomial, guesses)
        else:
            modes=zero_algo(self.characteristic_polynomial, f_starting_point, f_end_point)

        return modes
    def eigenvectors(self):
        """
        Computes the eigenvectors (null space) of the total admittance matrix
        at each eigenfrequency found in the range.

        Handles two cases:
        - Real eigenvalues (jω): eigenvalue is scalar
        - Complex eigenvalues (k + jω): eigenvalue is a tuple/list of length 2

        Returns:
            List of null space vectors (eigenvectors) for each mode.
        """
        #circuit_eigenvalues=self.eigenvalues(f_starting_point, f_end_point)
        circuit_eigenvalues=self.modes
        circuit_eigenvectors=[]

        if not circuit_eigenvalues:
            raise ValueError("No eigenvalues found.")
        
        #find max node in the circuit
        #max_node = max(max(comp.node_plus, comp.node_minus) for comp in self.netlist)
        components=self.netlist
        max_node = 0
        for comp in components:
            if hasattr(comp, "node_minus") and hasattr(comp, "node_plus"):
                max_node = max(max_node, comp.node_minus, comp.node_plus)
            elif hasattr(comp, "nodes"):
                max_node = max(max_node, *comp.nodes)
       # print(max_node)
        
        if isinstance(circuit_eigenvalues[0],numbers.Number):
            for eigen in circuit_eigenvalues:
                Y_f0=self.build_total_Y_matrix(1j*2*np.pi*1e9*eigen)
                
                if max_node==1:
                    full_vec=[0,1]
                else:

                    null_vecs = scipy.linalg.null_space(Y_f0,rcond=1e-10) #if it doesn't find the kernell, increase the tollerance rcond=1e-10
                    full_vec = np.zeros(max_node + 1, dtype=complex) #extent with ground value
                    full_vec[1:] = null_vecs.flatten()
                circuit_eigenvectors.append(full_vec)
        
        elif isinstance(circuit_eigenvalues[0], (list, tuple)) and len(circuit_eigenvalues[0]) == 2:
            for f, k in circuit_eigenvalues:
                z = 2 * np.pi * 1e6 * k + 1j * 2 * np.pi * 1e9 * f
                Y_z0=self.build_total_Y_matrix(z)
                if max_node==1:
                    full_vec=[0,1]
                else:
                    null_vecs = scipy.linalg.null_space(Y_z0,rcond=1e-10) #if it doesn't find the kernell, increase the tollerance
                    #da trattare il caso degenre??
                    full_vec = np.zeros(max_node + 1, dtype=complex) #extent with ground value
                    full_vec[1:] = null_vecs.flatten()
                circuit_eigenvectors.append(full_vec)

        else:
            raise ValueError("Eigenvalue format not recognized.")
        return circuit_eigenvectors
    
    def complex_frequencies(self):
        complex_f=[]
        #circuit_eigenvalues=self.eigenvalues(f_starting_point, f_end_point)
        circuit_eigenvalues=self.modes
        for eigen in circuit_eigenvalues:
            if isinstance(eigen,numbers.Number):
                complex_f.append(1j*2*np.pi*1e9*eigen)
            elif isinstance(eigen, (list, tuple)) and len(eigen) == 2:
                complex_f.append(2*np.pi*1e6*eigen[1]+1j*2*np.pi*1e9*eigen[0])
            else:
                raise ValueError("Eigenvalue format not recognized.")
        
        return complex_f

    def total_inductive_energy(self):
        """
        Calculate the total inductive energy stored into the circuit
        """
        circuit_eigenvalues=self.complex_frequencies()
        eigenvectors_with_ground=self.eigenvectors()
        E_inductive=[]

        for i in range(len(circuit_eigenvalues)):
            E_tot=0
            complex_f=circuit_eigenvalues[i]
            eigenvectors=eigenvectors_with_ground[i]
            for comp in self.netlist:
                if isinstance(comp, J):
                    current=comp.admittance(complex_f)*(eigenvectors[comp.node_plus]-eigenvectors[comp.node_minus])
                    E_tot+=comp.J_value*abs(current)**2/2
                if isinstance(comp, L):
                    current=comp.admittance(complex_f)*(eigenvectors[comp.node_plus]-eigenvectors[comp.node_minus])
                    E_tot+=comp.L_value*abs(current)**2/2
                if isinstance(comp,CPW):
                    #da verificare se è giusto !!!!!
                    E_tot+=comp.inductive_energy(eigenvectors[comp.node_plus],eigenvectors[comp.node_minus],complex_f)
                    #print('E cpw',comp.inductive_energy(eigenvectors[comp.node_plus],eigenvectors[comp.node_minus],complex_f))
                if isinstance(comp, CPW_coupler):
                    V=[eigenvectors[node] for node in comp.nodes]
                    E_tot+=comp.inductive_energy(V,complex_f)
                    #print('E coupler',comp.inductive_energy(V,complex_f))
            
            E_inductive.append(E_tot)
        return E_inductive



    def run_epr(self):
        """
        Implement energy participation ratio method to calculate Cross-Kerr matrix
        Returns:
            Cross-Kerr matrix
        """
        circuit_eigenvalues=self.complex_frequencies()
        f=[z.imag/2/np.pi for z in circuit_eigenvalues]
        eigenvectors_with_ground=self.eigenvectors()
        comp=self.netlist
        N_junct=0 #number of junction in the netlist
        junction_index=[] #index of the junction elements in the netlist

        #find junction
        for i in range(len(self.netlist)):
            if isinstance(self.netlist[i], J):
                N_junct+=1
                junction_index.append(i)
        
        if N_junct==0:
            raise ValueError('No junctions in the circuit')

        p=np.zeros((len(circuit_eigenvalues),N_junct)) #energy participation coefficients matrix
        E_tot=self.total_inductive_energy()

        #calculate energy participatio ratio matrix
        for m in range(len(circuit_eigenvalues)):
            for j in range(N_junct):
                eigenvectors=eigenvectors_with_ground[m]
                current=comp[junction_index[j]].admittance(circuit_eigenvalues[m])*(eigenvectors[comp[junction_index[j]].node_plus]-eigenvectors[comp[junction_index[j]].node_minus])
                Ej=comp[junction_index[j]].J_value*abs(current)**2/2
                p_mj=Ej/E_tot[m]
                p[m,j]=p_mj

        #calculate cross-kerr and self-kerr in matrix form
        chi=np.zeros((len(circuit_eigenvalues), len(circuit_eigenvalues)))
        for m in range(len(circuit_eigenvalues)):
            for n in range(len(circuit_eigenvalues)):
                for j in range(N_junct):
                    chi[m,n]+=0.25*f[m]*f[n]*p[m,j]*p[n,j]/comp[junction_index[j]].Ej()/comp[junction_index[j]].N**2/1e6  #MHZ 
        for m in range(len(circuit_eigenvalues)):
            chi[m,m]=chi[m,m]/2
        
        return chi

    def mode_frequencies(self):
        """
        Returns the frequencies of the modes of the circuit
        """
        #eigen=self.eigenvalues(f_starting_point, f_end_point)
        eigen=self.modes
        frequencies=[]

        if isinstance(eigen[0],numbers.Number):
            for val in eigen:
                frequencies.append(val)
        elif isinstance(eigen[0], (list, tuple)) and len(eigen[0]) == 2:
            for val in eigen:
                frequencies.append(val[0])
        return frequencies

    def kappa(self):
        """
        Returns the kappa of the modes of the circuit
        """
        #eigen=self.eigenvalues(f_starting_point, f_end_point)
        eigen=self.modes
        kappa=[]

        if isinstance(eigen[0],numbers.Number):
            for val in eigen:
                kappa.append(0)
        elif isinstance(eigen[0], (list, tuple)) and len(eigen[0]) == 2:
            for val in eigen:
                kappa.append(-2*val[1])
        return kappa
       
    def show_modes(self):
        """
        Function to visualize the modes of the circuit
        """
        #eigen=self.eigenvalues(f_starting_point, f_end_point)
        eigen=self.modes
        table=[]

        if isinstance(eigen[0],numbers.Number):
            for i, val in enumerate(eigen, 1):
                freq=val #GHz
                k=0 #Mhz
                table.append([i, f"{freq:.2e}", f"{k:.2e}"])
            print(tabulate(table, headers=["Mode", "Freq [GHz]", "k [MHz]"], tablefmt="pretty"))

        elif isinstance(eigen[0], (list, tuple)) and len(eigen[0]) == 2:
            for i, val in enumerate(eigen, 1):
                freq=val[0] #GHz
                k=2*val[1] #Mhz
                Q=freq*1e9/k/1e6
                table.append([i, f"{freq:.2e}", f"{k:.2e}", f"{Q:.2e}"])
            print(tabulate(table, headers=["Mode", "Freq [GHz]", "k [MHz]","Q"], tablefmt="pretty"))
            return
        
    def show_chi(self):
        """
        Function to visualize the Cross-Kerr matrix
        """
        chi=self.run_epr()
        N = chi.shape[0]
        table = []

        headers = ["Mode"] + [f"{j+1}" for j in range(N)]
        for i in range(N):
            row = [f"{i+1}"]
            for j in range(N):
                row.append(f"{chi[i,j]:.2e}")
            table.append(row)

        print("Chi matrix [MHz]:")
        print(tabulate(table, headers=headers, tablefmt="pretty"))
        return

    def show_all(self):
        """
        Show all the key parameter of the circuit
        """
        self.show_modes()
        self.show_chi()
        return

    def get_Z_submatrix(self,port,f,k=0):
        z=2*np.pi*1e6*k + 1j*2*np.pi*1e9*f
        Y=self.build_total_Y_matrix(z)
        Z=np.linalg.inv(Y)
        Z_submatrix=np.zeros((len(port),len(port)),dtype=complex)
        for i in range(len(port)):
            for j in range(len(port)):
                Z_submatrix[i,j]=Z[port[i]-1,port[j]-1]
        return Z_submatrix


