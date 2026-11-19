import numpy as np
from qutip import destroy, qeye, tensor
from copy import deepcopy
import numbers

import qultra.constants as constants
from ..core import QCircuit
from ..components.junction import J
from ..components.inductor import L
from ..components.cpw_transmission_line import CPW


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _has_flux(phi_ext):
    """True se phi_ext (numero o lista) contiene almeno un valore diverso da zero."""
    if phi_ext is None:
        return False
    return bool(np.any(np.asarray(phi_ext, dtype=float) != 0))


def _phase_fixed_eigenvectors(circuit: QCircuit):
    """
    Autovettori delle tensioni di nodo (con il ground) con la fase globale fissata:
    la componente di modulo massimo viene resa reale e positiva.
    Tutti i segni (giunzioni, induttori, CPW) sono calcolati da questi vettori,
    cosi' sono consistenti tra loro all'interno di ogni modo.
    """
    vecs = []
    for v in circuit.eigenvectors():
        v = np.asarray(v, dtype=complex)
        v = v * np.exp(-1j * np.angle(v[np.argmax(np.abs(v))]))
        vecs.append(v)
    return vecs


def _find_indices(circuit: QCircuit):
    """Indici nella netlist di: tutte le giunzioni, induttori con flusso, CPW con flusso."""
    junction_index, inductor_index, cpw_index = [], [], []
    for i, comp in enumerate(circuit.netlist):
        if isinstance(comp, J):
            junction_index.append(i)
        elif isinstance(comp, L) and _has_flux(comp.phi_ext):
            inductor_index.append(i)
        elif isinstance(comp, CPW) and _has_flux(comp.phi_ext):
            if len(comp.phi_ext) != 2:
                raise ValueError("phi_ext of a CPW must be [phase at x=0 (node_plus), phase at x=l (node_minus)]")
            cpw_index.append(i)
    return junction_index, inductor_index, cpw_index


# ---------------------------------------------------------------------------
# Matrici dei segni (lumped)
# ---------------------------------------------------------------------------

def _sign_matrix(circuit: QCircuit, index, vecs):
    """
    Segno del modo m sul ramo index[j]: sign(Re(V_plus - V_minus)) con fase fissata.
    Per un ramo induttivo V = j*omega*L*I, quindi e' anche il segno del flusso di ramo.
    """
    comp = circuit.netlist
    s = np.zeros((len(vecs), len(index)))
    for m, v in enumerate(vecs):
        for j, k in enumerate(index):
            V = v[comp[k].node_plus] - v[comp[k].node_minus]
            s[m, j] = 1.0 if V.real >= 0 else -1.0   # evita sign(0)=0
    return s


# ---------------------------------------------------------------------------
# EPR per induttori con flusso
# ---------------------------------------------------------------------------

def _inductor_epr(circuit: QCircuit, inductor_index, vecs, E_tot):
    """
    Energy participation ratio degli induttori con flusso esterno.
    Convenzione: E = 1/2 L |I|^2 (energia di picco), la stessa di total_inductive_energy().
    """
    lam = circuit.complex_frequencies()
    comp = circuit.netlist
    p = np.zeros((len(lam), len(inductor_index)))
    for m, z in enumerate(lam):
        v = vecs[m]
        for j, k in enumerate(inductor_index):
            current = comp[k].admittance(z) * (v[comp[k].node_plus] - v[comp[k].node_minus])
            p[m, j] = comp[k].L_value * abs(current)**2 / 2 / E_tot[m]
    return p


# ---------------------------------------------------------------------------
# Correnti di punto zero agli estremi delle CPW con flusso
# ---------------------------------------------------------------------------

def _cpw_zpf_currents(circuit: QCircuit, cpw_index, vecs, E_tot):
    """
    Correnti di punto zero con segno [A] agli estremi delle CPW.

    I0[m, j] : corrente in x = 0 (node_plus)
    Il[m, j] : corrente in x = l (node_minus)
    Verso positivo: da node_plus verso node_minus (convenzione di CPW.current).

    I_zpf = |I| * sqrt(hbar*omega / (4*E_tot)), con E_tot energia induttiva di picco
    (1/2 L|I|^2, CPW incluse). Il segno e' sign(Re(j*I)), coerente con sign(Re(V))
    usato per giunzioni e induttori (per un induttore j*omega*L*I = V).
    """
    lam = circuit.complex_frequencies()
    comp = circuit.netlist
    I0 = np.zeros((len(lam), len(cpw_index)))
    Il = np.zeros((len(lam), len(cpw_index)))
    for m, z in enumerate(lam):
        v = vecs[m]
        norm = np.sqrt(constants.hbar * z.imag / (4 * E_tot[m]))
        for j, k in enumerate(cpw_index):
            c = comp[k]
            for arr, x in ((I0, 0.0), (Il, c.l)):
                I = c.current(v[c.node_plus], v[c.node_minus], z, x)
                sign = 1.0 if (1j * I).real >= 0 else -1.0
                arr[m, j] = norm * abs(I) * sign
    return I0, Il


# ---------------------------------------------------------------------------
# Hamiltoniana
# ---------------------------------------------------------------------------

def fluxes_hamiltonian(circuit: QCircuit, excitations):
    """
    Hamiltoniana [Hz] del circuito con flussi esterni.

    Convenzioni sui flussi esterni (fasi adimensionali 2*pi*Phi/Phi0, gia' moltiplicate per b):
      - L.phi_ext   : numero
      - J.phi_ext   : numero (entra come cos(phi/N + phi_ext))
      - CPW.phi_ext : [fase all'estremo x=0 (node_plus), fase all'estremo x=l (node_minus)]
    """
    modes = circuit.mode_frequencies()  # GHz

    if isinstance(excitations, numbers.Number):
        number_of_excitations = [excitations for _ in range(len(modes))]
    elif isinstance(excitations, list):
        if len(excitations) != len(modes):
            raise ValueError("The length of the list of excitations must be equal to the number of modes")
        number_of_excitations = excitations
    else:
        raise ValueError("excitations must be a number or a list of numbers")

    comp = circuit.netlist
    junction_index, inductor_index, cpw_index = _find_indices(circuit)

    # quantita' calcolate una sola volta
    vecs = _phase_fixed_eigenvectors(circuit)
    E_tot = circuit.total_inductive_energy()

    scaling = (constants.phi0 / (2 * np.pi))**2 / constants.h   # [J/H -> Hz]: scaling/L = E_L in Hz
    scaling_cpw = (constants.phi0 / (2 * np.pi)) / constants.h  # Phi0/2pi * I [A] -> Hz

    # operatori di distruzione
    qeye_list = [qeye(n) for n in number_of_excitations]
    operators = []
    H_lin = 0
    for index, mode in enumerate(modes):
        a_to_tensor = deepcopy(qeye_list)
        a_to_tensor[index] = destroy(number_of_excitations[index])
        a = tensor(a_to_tensor)
        operators.append(a)
        H_lin += 1e9 * mode * a.dag() * a  # Hz

    # --- induttori con flusso: E_L * phi_ext * phi_L
    if inductor_index:
        p_inductor = _inductor_epr(circuit, inductor_index, vecs, E_tot)
        s_inductor = _sign_matrix(circuit, inductor_index, vecs)
        for j, k in enumerate(inductor_index):
            phi_inductor = 0
            for m, a in enumerate(operators):
                phi_inductor += s_inductor[m, j] * np.sqrt(
                    p_inductor[m, j] * modes[m] * 1e9 / (2 * comp[k].El())) * (a + a.dag())
            H_lin += (scaling / comp[k].L_value) * comp[k].phi_ext * phi_inductor

    # --- CPW con flusso: (Phi0/2pi) * phi_ext * I_zpf * (a + a^dag) a ciascun estremo
    if cpw_index:
        I0, Il = _cpw_zpf_currents(circuit, cpw_index, vecs, E_tot)
        for j, k in enumerate(cpw_index):
            phi_0, phi_l = comp[k].phi_ext
            for m, a in enumerate(operators):
                H_lin += scaling_cpw * (phi_0 * I0[m, j] + phi_l * Il[m, j]) * (a + a.dag())

    if not junction_index:
        return 0.5 * (H_lin + H_lin.dag())  # elimina residui numerici non hermitiani

    # --- giunzioni
    _, p = circuit.run_epr()
    s = _sign_matrix(circuit, junction_index, vecs)

    phi = [0 * operators[0] for _ in junction_index]
    for j, k in enumerate(junction_index):
        for m, a in enumerate(operators):
            phi[j] += s[m, j] * np.sqrt(p[m, j] * modes[m] * 1e9 / (2 * comp[k].Ej())) * (a + a.dag())

    # parte non lineare: N^2 Ej cos(phi/N + phi_ext) + Ej phi^2/2 (sottrae il quadratico gia' in H_lin)
    H_nl = 0
    for j, k in enumerate(junction_index):
        Jc = comp[k]
        H_nl += Jc.N**2 * Jc.Ej() * (np.exp(1j * Jc.phi_ext) * (1j * phi[j] / Jc.N).expm()
                                     + np.exp(-1j * Jc.phi_ext) * (-1j * phi[j] / Jc.N).expm()) / 2
        H_nl += Jc.Ej() * phi[j]**2 / 2

    H_total = H_lin - H_nl
    return 0.5 * (H_total + H_total.dag())  # elimina residui numerici non hermitiani