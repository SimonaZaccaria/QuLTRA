import numpy as np
from qutip import destroy, qeye, tensor
from copy import deepcopy
import numbers

import qultra.constants as constants
from ..core import QCircuit
from ..components.junction import J
from ..components.inductor import L


def _sign_matrix_junction(circuit: QCircuit, p):
    '''
    sm=-1*np.ones((p.shape[0],p.shape[1]))
    for m in range(p.shape[0]):
        sm[m,np.argmax(p[m,:])]=1 #creo la matrice dei segni
    for j in range(p.shape[1]):
        sm[np.argmax(p[:,j]),j]=1 #creo la matrice dei segni#creo la matrice dei segni
    '''
    circuit_eigenvalues=circuit.complex_frequencies()
    f=[z.imag/2/np.pi for z in circuit_eigenvalues]
    eigenvectors_with_ground=circuit.eigenvectors()
    comp=circuit.netlist
    N_junct=0 #number of junction in the netlist
    junction_index=[] #index of the junction elements in the netlist

    #find junction
    for i in range(len(circuit.netlist)):
        if isinstance(circuit.netlist[i], J):
            N_junct+=1
            junction_index.append(i)

    if N_junct==0:
        raise ValueError('No junctions in the circuit')

    s=np.zeros((len(circuit_eigenvalues),N_junct))
    #calculate energy participatio ratio matrix
    for m in range(len(circuit_eigenvalues)):
        for j in range(N_junct):
            eigenvectors=eigenvectors_with_ground[m]
            V=(eigenvectors[comp[junction_index[j]].node_plus]-eigenvectors[comp[junction_index[j]].node_minus])
            s[m,j]=np.sign(V.real)
    #print(s)
    return s

def _sign_matrix_inductor(circuit: QCircuit, p):
    '''
    sm=-1*np.ones((p.shape[0],p.shape[1]))
    for m in range(p.shape[0]):
        sm[m,np.argmax(p[m,:])]=1 #creo la matrice dei segni
    for j in range(p.shape[1]):
        sm[np.argmax(p[:,j]),j]=1 #creo la matrice dei segni#creo la matrice dei segni
    '''
    circuit_eigenvalues=circuit.complex_frequencies()
    f=[z.imag/2/np.pi for z in circuit_eigenvalues]
    eigenvectors_with_ground=circuit.eigenvectors()
    comp=circuit.netlist
    N_inductor=0 #number of inductor in the netlist
    inductor_index=[] #index of the inductor elements in the netlist

    #find inductor
    for i in range(len(circuit.netlist)):
        if isinstance(circuit.netlist[i], L):
            N_inductor+=1
            inductor_index.append(i)

    if N_inductor==0:
        raise ValueError('No inductors in the circuit')

    s=np.zeros((len(circuit_eigenvalues),N_inductor))
    #calculate energy participatio ratio matrix
    for m in range(len(circuit_eigenvalues)):
        for j in range(N_inductor):
            eigenvectors=eigenvectors_with_ground[m]
            V=(eigenvectors[comp[inductor_index[j]].node_plus]-eigenvectors[comp[inductor_index[j]].node_minus])
            s[m,j]=np.sign(V.real)
    #print(s)
    return s
def _inductor_epr(circuit: QCircuit):
    """
    Compute the Cross-Kerr matrix using the energy participation ratio method.

    Returns
    -------
    chi: numpy.ndarray
        The Cross-Kerr matrix of the system [MHz].
    p: numpy.ndarray
        The energy participation ratio matrix.
    """
    circuit_eigenvalues=circuit.complex_frequencies()
    f=[z.imag/2/np.pi for z in circuit_eigenvalues]
    eigenvectors_with_ground=circuit.eigenvectors()
    comp=circuit.netlist
    N_inductor=0 #number of inductor in the netlist
    inductor_index=[] #index of the inductor elements in the netlist

    #find inductor
    for i in range(len(circuit.netlist)):
        if isinstance(circuit.netlist[i], L):
            if circuit.netlist[i].phi_ext!=0:
                N_inductor+=1
                inductor_index.append(i)

    if N_inductor==0:
        raise ValueError('No inductors in the circuit')

    p=np.zeros((len(circuit_eigenvalues),N_inductor)) #energy participation coefficients matrix
    
    
    E_tot=circuit.total_inductive_energy()

    #calculate energy participatio ratio matrix
    for m in range(len(circuit_eigenvalues)):
        for j in range(N_inductor):
            eigenvectors=eigenvectors_with_ground[m]
            current=comp[inductor_index[j]].admittance(circuit_eigenvalues[m])*(eigenvectors[comp[inductor_index[j]].node_plus]-eigenvectors[comp[inductor_index[j]].node_minus])
            El=comp[inductor_index[j]].L_value*abs(current)**2/2
            p_mj=El/E_tot[m]
            p[m,j]=p_mj
    return p

def fluxes_hamiltonian(circuit: QCircuit,excitations):

    modes=circuit.mode_frequencies() #GHz
    

    if isinstance(excitations,numbers.Number):
        number_of_excitations=[excitations for _ in range(len(modes))]
    elif isinstance(excitations,list):
        if len(excitations)!=len(modes):
            raise ValueError("The length of the list of excitations must be equal to the number of modes")
        number_of_excitations=excitations
    else:
        raise ValueError("excitations must be a number or a list of numbers")
    
    
    comp=circuit.netlist
    N_junct=0 #number of junction in the netlist
    N_inductance=0 #number of inductance in the netlist
    junction_index=[]
    inductor_index=[]
    qeye_list = [qeye(n) for n in number_of_excitations]
    operators=[]
    scaling=(constants.phi0/(2*np.pi))**2/constants.h #scaling factor to express the Hamiltonian in Hz

    #find junction
    for i in range(len(circuit.netlist)):
        if isinstance(circuit.netlist[i], J):
            N_junct+=1
            junction_index.append(i)


    #find inductor
    for i in range(len(circuit.netlist)):
        if isinstance(circuit.netlist[i], L):
            if circuit.netlist[i].phi_ext!=0:
                N_inductance+=1
                inductor_index.append(i)

    #create the linear part of the Hamiltonian
    #H_lin = tensor([qeye(n) for n in excitations]) * 0
    H_lin=0
    for index,mode in enumerate(modes):
        a_to_tensor = deepcopy(qeye_list)
        a_to_tensor[index] = destroy(number_of_excitations[index])
        a = tensor(a_to_tensor)
        operators.append(a)
        H_lin += 1e9*mode*a.dag()*a #Hamiltonian expressed in Hz
    
    if N_inductance!=0:
        p_inductor=_inductor_epr(circuit)
        s_inductor=_sign_matrix_inductor(circuit,p_inductor)
        
        #phi_inductor = [0 * operators[0] for _ in range(N_inductance)]
        for j in range(N_inductance):
            phi_inductor=0
            for m in range(len(modes)):
                a = operators[m]
                phi_inductor+=s_inductor[m,j]*np.sqrt(p_inductor[m,j]*modes[m]*1e9/(2*comp[inductor_index[j]].El()))* (a + a.dag()) #zpf by using epr formula
            H_lin+=(scaling/comp[inductor_index[j]].L_value)*comp[inductor_index[j]].phi_ext*phi_inductor #add the external flux contribution to the linear part of the Hamiltonian

    if N_junct==0:
        return H_lin
    
    _,p=circuit.run_epr() #participation ratio coefficients
    s=_sign_matrix_junction(circuit,p) #sign matrix
    #calculate phi operators
    phi = [0 * operators[0] for _ in range(N_junct)]

    '''
    for m in range(len(modes)):
        a = operators[m]
        for j in range(N_junct):
            phi[j]+=np.sqrt(p[m,j]*modes[m]*1e9/(2*comp[junction_index[j]].Ej()))* (a + a.dag()) #zpf by using epr formula
    '''
    for j in range(N_junct):
        for m in range(len(modes)):
            a = operators[m]
            phi[j]+=s[m,j]*np.sqrt(p[m,j]*modes[m]*1e9/(2*comp[junction_index[j]].Ej()))* (a + a.dag()) #zpf by using epr formula
            
    
    #create nonlinear part of the Hamiltonian
    #H_nl = tensor([qeye(n) for n in excitations]) *0
    H_nl=0
    

    for j in range(N_junct):
        #da chiedere il verso del flusso esterno
        
        H_nl+=(comp[junction_index[j]].N**2* comp[junction_index[j]].Ej()*(np.exp(1j*comp[junction_index[j]].phi_ext)*(1j*phi[j]/comp[junction_index[j]].N).expm()+np.exp(-1j*comp[junction_index[j]].phi_ext)*(-1j*phi[j]/comp[junction_index[j]].N).expm())/2)+comp[junction_index[j]].Ej()*phi[j]**2/2

        #H_nl+=comp[junction_index[j]].N**2* comp[junction_index[j]].Ej()*(phi[j]/comp[junction_index[j]].N).cosm()+ comp[junction_index[j]].Ej()*phi[j]**2/2
    
    H_total=H_lin-H_nl
    
    return H_total