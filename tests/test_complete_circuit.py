import numpy as np
import qultra as qu

def test_complete_circuit():
    Cj=90e-15 #qubit capacitance
    Lj=8e-9 #qubit junction
    Cg=5e-15 #coupling capacitance
    l=4.5e-3 #cpw length

    l_coupler=0.7e-3 #coupler length
    coupler_nodes=[3,0,4,5] #coupler nodes
    gap=[10,10,10,10] #coupler gaps
    width=[15,10,15] #coupler widths

    f_min=3 #minimum frequency [GHz]
    f_max=12 #maximum frequency [GHz]

    net1=[qu.C(0,1,Cj),qu.J(0,1,Lj),qu.C(1,2,Cg),qu.CPW(3,2,l-l_coupler),qu.CPW_coupler(coupler_nodes,gap,width,l_coupler),qu.R(0,4,50),qu.R(0,5,50)]

    circuit_with_loss=qu.QCircuit(net1,f_min,f_max)

    f= circuit_with_loss.mode_frequencies()
    k= circuit_with_loss.kappa()
    chi,_= circuit_with_loss.run_epr()

    assert isinstance(f, list)
    assert isinstance(k, list)
    assert isinstance(chi, np.ndarray)
    assert len(f)==2
    assert chi.shape[0]==chi.shape[1]

