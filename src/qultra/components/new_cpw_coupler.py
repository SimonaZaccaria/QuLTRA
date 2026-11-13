import numpy as np
import qultra.constants as constants
from ..cm_flipchip import G_matrices, CL_matrices
from .cpw_coupler import CPW_coupler


class new_CPW_coupler(CPW_coupler):
    """
    CPW coupler (coupled lines) also in flip-chip architecture.
    The propagation is described with the modes of the coupled lines,
    each one with its own phase velocity.

    Node ordering (top view)
     ::

        1----------2

        3----------4

    Parameters
    ----------
    nodes : list of int
        List of 4 nodes [start line 1, end line 1, start line 2, end line 2].
    gap : list of float
        List of gaps [um], len(gap) = len(cpw) + 1.
    cpw : list of float
        List of strip widths [um].
    l : float
        Length of the coupler [m].
    hs : float, optional
        Distance between the two chips in flip-chip [um]. Default is None (single chip).

    Notes
    -----
    l is in meters, while gap, cpw and hs are in micrometers.
    constants.epsilon_r is read when the component is created.
    """

    def __init__(self, nodes, gap, cpw, l, hs=None):
        if len(nodes) != 4:
            raise ValueError("The component must be connected to 4 nodes")
        if len(gap) != len(cpw) + 1:
            raise ValueError(
                f"`gap` must have length equal to `len(cpw) + 1`. Got len(gap) = {len(gap)}, len(cpw) = {len(cpw)}"
            )
        self.nodes = nodes
        self.gap = gap
        self.cpw = cpw
        self.l = l        # [m]
        self.hs = hs      # [um]

        # matrici per unita' di lunghezza
        G_top, G_bottom = G_matrices(gap, cpw, hs)
        self.C, self.L = CL_matrices(G_top, G_bottom)

        # modi: autovettori T (colonne) e autovalori 1/v_m^2 di L C
        lam, T = np.linalg.eig(self.L @ self.C)
        self.T = T.real
        self.Tinv = np.linalg.inv(self.T)
        self.v_modes = 1 / np.sqrt(lam.real)
        self.epsilon_eff_modes = (constants.c / self.v_modes)**2

        # ammettenza caratteristica in forma modale: ampiezze modali -> correnti delle strisce
        self.Yc_modal = self.C @ self.T @ np.diag(self.v_modes)

    def _modal_amplitudes(self, Vm_0, Vm_l, z):
        """Ampiezze progressive e regressive di ogni modo, date le tensioni modali agli estremi."""
        exp_z = np.exp(self.l * z / self.v_modes)
        exp_mz = np.exp(-self.l * z / self.v_modes)
        imag_part = exp_z - exp_mz
        V_plus = (Vm_0 * exp_z - Vm_l) / imag_part
        V_minus = (-Vm_0 * exp_mz + Vm_l) / imag_part
        return V_plus, V_minus, exp_z, exp_mz

    def Y(self, z):
        """
        Construct the complex admittance matrix of the coupler.

        Parameters
        ----------
        z : complex
            Complex variable (complex frequency)

        Returns
        -------
        Y_matrix : numpy array
            complex admittance matrix (node ordering of self.nodes)
        """
        dim = self.C.shape[0]
        Y_matrix = np.zeros((2 * dim, 2 * dim), dtype=np.complex128)
        zero = np.zeros(dim)
        for i in range(2 * dim):
            k = i // 2
            Vm = self.Tinv[:, k]                   # 1 V sulla linea k, scomposto nei modi
            if i % 2 == 0:                         # 1 V all'inizio della linea k
                V_plus, V_minus, exp_z, exp_mz = self._modal_amplitudes(Vm, zero, z)
            else:                                  # 1 V alla fine della linea k
                V_plus, V_minus, exp_z, exp_mz = self._modal_amplitudes(zero, Vm, z)
            I_start = self.Yc_modal @ (V_plus - V_minus)
            I_end = self.Yc_modal @ (V_plus * exp_mz - V_minus * exp_z)
            for j in range(dim):
                Y_matrix[2 * j, i] = I_start[j]
                Y_matrix[2 * j + 1, i] = -I_end[j]
        return Y_matrix

    def current(self, V, z, x):
        """
        Calculate the currents in the coupler lines given the node voltages.

        Parameters
        ----------
        V : complex array
            Voltages at the nodes (same ordering as self.nodes).
        z : complex
            Complex variable (complex frequency).
        x : float
            Position along the coupler [m].

        Returns
        -------
        I : complex array, shape (dim, 1)
            Currents in the lines at position x.
        """
        V = np.asarray(V, dtype=np.complex128).ravel()
        Vm_0 = self.Tinv @ V[0::2]                 # tensioni modali all'inizio
        Vm_l = self.Tinv @ V[1::2]                 # tensioni modali alla fine
        V_plus, V_minus, _, _ = self._modal_amplitudes(Vm_0, Vm_l, z)
        I = self.Yc_modal @ (V_plus * np.exp(-x * z / self.v_modes) - V_minus * np.exp(x * z / self.v_modes))
        return I.reshape(-1, 1)

    # inductive_energy e' ereditato da CPW_coupler: integra I^H L I / 2 lungo il tratto,
    # usando self.current e self.L, che qui sono gia' quelli corretti.