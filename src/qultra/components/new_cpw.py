import qultra.constants as constants
from ..cm_flipchip import cpw_parameters
from .cpw_transmission_line import CPW


class new_CPW(CPW):
    """
    CPW defined by its geometry, also in flip-chip architecture.
    Characteristic impedance and phase velocity are computed from the
    cross section with conformal mapping.

    Parameters
    ----------
    node_minus : int
        The node to which the negative terminal is connected
    node_plus : int
        The node to which the positive terminal is connected
    l : float
        Length of the CPW [m]
    w : float
        Width of the center conductor [um]
    g : float
        Gap between the center conductor and the ground planes [um]
    hs : float, optional
        Distance between the two chips in flip-chip [um]. Default is None (single chip).
    phi_ext : list of float, optional
        External applied phase at the two ends of the line. Default is None.

    Notes
    -----
    l is in meters, while w, g and hs are in micrometers.
    constants.epsilon_r is read when the component is created.
    """
    def __init__(self, node_minus, node_plus, l, w, g, hs=None, phi_ext=None):
        self.w = w      # [um]
        self.g = g      # [um]
        self.hs = hs    # [um]
        Z0, v = cpw_parameters(w, g, hs)
        super().__init__(node_minus, node_plus, l, Z0=Z0, phi_ext=phi_ext)
        self.v = v                                   # sovrascrive la v calcolata da CPW
        self.epsilon_eff = (constants.c / v)**2      # coerente con l'attributo di CPW