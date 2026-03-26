import numpy as np
from qutip import destroy, qeye, tensor
from copy import deepcopy
import numbers
from math import factorial

from ..core import QCircuit
from ..components.junction import J


def _sign_matrix(circuit: QCircuit, p):
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


def nofluxes_hamiltonian(circuit: QCircuit,excitations,taylor=False,order=4):

    modes=circuit.mode_frequencies() #GHz
    

    if isinstance(excitations,numbers.Number):
        number_of_excitations=[excitations for _ in range(len(modes))]
    elif isinstance(excitations,list):
        if len(excitations)!=len(modes):
            raise ValueError("The length of the list of excitations must be equal to the number of modes")
        number_of_excitations=excitations
    else:
        raise ValueError("excitations must be a number or a list of numbers")
    
    if order<4:
        raise ValueError("The order of the Taylor expansion must be >=4")
    
    comp=circuit.netlist
    N_junct=0 #number of junction in the netlist
    junction_index=[]
    qeye_list = [qeye(n) for n in number_of_excitations]
    operators=[]

    #find junction
    for i in range(len(circuit.netlist)):
        if isinstance(circuit.netlist[i], J):
            N_junct+=1
            junction_index.append(i)

    
    #create the linear part of the Hamiltonian
    #H_lin = tensor([qeye(n) for n in excitations]) * 0
    H_lin=0
    for index,mode in enumerate(modes):
        a_to_tensor = deepcopy(qeye_list)
        a_to_tensor[index] = destroy(number_of_excitations[index])
        a = tensor(a_to_tensor)
        operators.append(a)
        H_lin += 1e9*mode*a.dag()*a #Hamiltonian expressed in Hz
    if N_junct==0:
        return H_lin
    
    _,p=circuit.run_epr() #participation ratio coefficients
    s=_sign_matrix(circuit,p) #sign matrix
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
    
    if taylor:
        #approximate the cosine potential with Taylor expansion
        for j in range(N_junct):
            expansion=0*operators[0]
            for n in range(2,order//2+1):
                expansion+=((-1)**n/factorial(2*n))*(phi[j]/comp[junction_index[j]].N)**(2*n)
            H_nl+=comp[junction_index[j]].N**2* comp[junction_index[j]].Ej()*expansion
    else:
        #use the exact cosine potential by expressing it thorugh Euler's exponentials
        for j in range(N_junct):
            #H_nl+=comp[junction_index[j]].N**2* comp[junction_index[j]].Ej()*((1j*phi[j]/comp[junction_index[j]].N).expm()+(-1j*phi[j]/comp[junction_index[j]].N).expm())/2
            H_nl+=comp[junction_index[j]].N**2* comp[junction_index[j]].Ej()*(phi[j]/comp[junction_index[j]].N).cosm()+ comp[junction_index[j]].Ej()*phi[j]**2/2
        
    H_total=H_lin-H_nl
    
    return H_total