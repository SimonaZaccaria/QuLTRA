import numpy as np
import scipy
import scipy.integrate
from scipy.optimize import root
from scipy.constants import epsilon_0

try:
    from . import constants          # import relativo, come negli altri file di QuLTRA
except ImportError:
    import constants                 # esecuzione dal sorgente senza installazione


def branch_point_coordinates(gap, cpw):
        ''' '''
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

def exponential_coordinates(a_coordinates,b_coordinates,shift):
        x0 = b_coordinates[-1].real / 2                     # origine al centro della struttura (precisione numerica) non utile concettualmente parlando
        a_exp = [complex(np.exp(np.pi*(x.real - x0)/shift)) for x in a_coordinates]
        b_exp = [complex(np.exp(np.pi*(x.real - x0)/shift)) for x in b_coordinates]
        return a_exp,b_exp

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
    if len(c_coordinates_guest)==0:    # una sola striscia: nessun gap tra masse,
        return []                      # quindi nessun punto c da cercare
    sol=root(c_solve,c_coordinates_guest)
    if not sol.success:
        raise RuntimeError("Root finding failed: " + sol.message)
    c_coordinates=[complex(x) for x in sol.x]
    return c_coordinates          

def G_matrices(gap, cpw, hs=None):
    """Fattori geometrici delle due meta' della sezione.
    G_top:    cavita' d'aria alta hs (hs=None: semispazio, uguale a G_bottom)
    G_bottom: semispazio inferiore (substrato)
    gap, cpw, hs nella stessa unita' (micrometri)."""
    a, b = branch_point_coordinates(gap, cpw)
    if hs is not None:
        a_exp, b_exp = exponential_coordinates(a, b, hs)
    n = len(cpw)
    G_bottom = np.zeros((n, n))
    G_top = np.zeros((n, n))
    for j in range(n):                                   # striscia j a 1 V
        ap, bp = conformal_mapping(a, b, find_c(a, b, j))
        if hs is not None:
            ap_exp, bp_exp = conformal_mapping(a_exp, b_exp, find_c(a_exp, b_exp, j))
        for i in range(n):
            G_bottom[i, j] = (bp[i].real - ap[i+1].real) / bp[j].imag
            if hs is None:
                G_top[i, j] = G_bottom[i, j]
            else:
                G_top[i, j] = (bp_exp[i].real - ap_exp[i+1].real) / bp_exp[j].imag
    if n == 3:                                           # massa centrale: messa a 0 V
        G_bottom = np.delete(np.delete(G_bottom, 1, axis=0), 1, axis=1)
        G_top = np.delete(np.delete(G_top, 1, axis=0), 1, axis=1)
    return G_top, G_bottom

def CL_matrices(G_top, G_bottom, epsilon_r=None):
    """Matrici di capacita' C' e induttanza L' per unita' di lunghezza [F/m], [H/m].
    epsilon_r=None: usa constants.epsilon_r, letto al momento della chiamata."""
    if epsilon_r is None:
        epsilon_r = constants.epsilon_r
    C_matrix = epsilon_0 * (epsilon_r * G_bottom + G_top)
    Ca_matrix = epsilon_0 * (G_bottom + G_top)
    L_matrix = np.linalg.inv(Ca_matrix) / constants.c**2
    return C_matrix, L_matrix

def cpw_parameters(w, g, hs=None):
    """Impedenza caratteristica [Ohm] e velocita' di fase [m/s] di una CPW.
    w, g, hs in micrometri; hs=None: chip singolo."""
    G_top, G_bottom = G_matrices([g, g], [w], hs)
    C_matrix, L_matrix = CL_matrices(G_top, G_bottom)
    C = C_matrix[0, 0]
    L = L_matrix[0, 0]
    Z0 = np.sqrt(L / C)
    v = 1 / np.sqrt(L * C)
    return Z0, v