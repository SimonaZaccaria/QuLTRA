import numpy as np
import scipy
import numbers
from tabulate import tabulate


from .components.capacitor import C 
from .components.inductor import L
from .components.resistor import R 
from .components.junction import J
from .components.cpw_transmission_line import CPW
from .components.cpw_coupler import CPW_coupler




try:
    from .find_zeros import*
except ImportError:
    # When running from source without pip installation
    from find_zeros import *




  

class QCircuit:
    """
    Represents a quantum circuit.

    Parameters
    ----------
    netlist : list
        List of component instances that define the circuit.
    f_starting_point : float
        The start frequency of the interval over which the circuit is analyzed.
    f_end_point : float
        The end frequency of the interval over which the circuit is analyzed.

    Attributes
    ----------
    modes : list of lists of length 2
        Each element is a list of two values:
            - The first value represents the eigenmode frequency in GHz.
            - The second value represents the eigenmode dissipation rate in MHz.
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
        ''' '''
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
        ''' '''
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
        
        Parameters
        ----------
            z: complex 
                complex variable (e.g., z = jω or k + jω)
        
        Returns
        -------
            Y_reduced: numpy matrix
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
        ''' '''
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
        ''' '''
        Y_matrix=self.build_total_Y_matrix(z)
        det_Y=np.linalg.det(Y_matrix)
        K = scipy.linalg.null_space(Y_matrix,rcond=1e-10)
        dim_kernel = K.shape[1]
        return det_Y, Y_matrix.shape[0],dim_kernel
    
    def there_is_R(self):
        ''' '''
        """
        Check if there are resistive components
        """
        for comp in self.netlist:
            if isinstance(comp, R):
                return True
            
        return False
    
    def eigenvalues(self,f_starting_point, f_end_point):
        ''' '''
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
        ''' '''
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
        ''' '''
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
        ''' '''
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
        Compute the Cross-Kerr matrix using the energy participation ratio method.

        Returns
        -------
        chi: numpy.ndarray
            The Cross-Kerr matrix of the system [MHz].
        p: numpy.ndarray
            The energy participation ratio matrix.
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
        
        return chi,p

    def mode_frequencies(self):
        """
        Returns the frequencies of the modes of the circuit [GHz]

        Returns
        -------
        frequencies: list
            Mode frequencies in GHz.
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
        Returns the kappa of the modes of the circuit [MHz]

        Returns
        -------
        kappa: list
            Mode dissipation rates in MHz.
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
        chi, _=self.run_epr()
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
        """
        Compute the impedance submatrix corresponding to a given set of nodes.

        Parameters
        ----------
        port : list of int
            List of node indices (1-based numbering) for which the impedance 
            submatrix is extracted.
        f : float
            Frequency in GHz at which the impedance matrix is computed.
        k : float, optional
            Dissipation rate in MHz at which the impedance matrix is computed. Default is 0.

        Returns
        -------
        Z_submatrix : ndarray of complex, shape (len(nodes), len(nodes))
            Impedance submatrix corresponding to the selected nodes.
        """
        z=2*np.pi*1e6*k + 1j*2*np.pi*1e9*f
        Y=self.build_total_Y_matrix(z)
        Z=np.linalg.inv(Y)
        Z_submatrix=np.zeros((len(port),len(port)),dtype=complex)
        for i in range(len(port)):
            for j in range(len(port)):
                Z_submatrix[i,j]=Z[port[i]-1,port[j]-1]
        return Z_submatrix
    
   
    def _there_are_fluxes(self):
        ''' '''
        """
        Check if there are fluxes in the circuit
        """
        for comp in self.netlist:
            if isinstance(comp, J):
                if comp.phi_ext!=0:
                    return True
            if isinstance(comp,L):
                if comp.phi_ext!=0:
                     return True
        return False
            
         
    def hamiltonian(self,excitations,taylor=False,order=4):
        from .simulations.nofluxes_simulation import nofluxes_hamiltonian
        from .simulations.fluxes_simulation import fluxes_hamiltonian

        if self._there_are_fluxes():
            return fluxes_hamiltonian(self,excitations)
        else:  
            return nofluxes_hamiltonian(self,excitations,taylor,order)



