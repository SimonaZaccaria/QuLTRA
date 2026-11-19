import numpy as np
from copy import deepcopy
import networkx as nx
import scipy.linalg as sl

from ..core import QCircuit
from ..components.junction import J
from ..components.inductor import L
from ..components.capacitor import C
from ..components.cpw_transmission_line import CPW
from ..components.cpw_coupler import CPW_coupler
from ..components.resistor import R


# ---------------------------------------------------------------------------
# Netlist equivalente
# ---------------------------------------------------------------------------

def _there_is_R(netlist):
    """
    Check if there are resistive components
    """
    for comp in netlist:
        if isinstance(comp, R):
            return True
    return False


def _reduce_netlist(netlist):
    """
    Remove resistive components from the netlist
    """
    return [comp for comp in netlist if not isinstance(comp, R)]


def _max_node(netlist):
    """
    Compute the max node of the network. The max node is the highest node number in the netlist.
    """
    max_node = 0
    for comp in netlist:
        if hasattr(comp, "node_minus") and hasattr(comp, "node_plus"):
            max_node = max(max_node, comp.node_minus, comp.node_plus)
        elif hasattr(comp, "nodes"):
            max_node = max(max_node, *comp.nodes)
    return max_node


def _equivalent_CPW_circuit(comp, new_node):
    """
    T equivalente di una CPW: L(node_minus, new) - C(new, 0) - L(new, node_plus).
    I valori sono quelli fisici della linea (i b non dipendono da essi, ma cosi' il sistema e' ben condizionato).
    Il primo induttore sta all'estremo x=l (node_minus), il secondo all'estremo x=0 (node_plus).
    """
    L_half = comp.Z0 * comp.l / comp.v / 2
    C_line = comp.l / (comp.Z0 * comp.v)
    return [L(comp.node_minus, new_node, L_half), C(new_node, 0, C_line), L(new_node, comp.node_plus, L_half)]


def _equivalent_CPW_coupler_circuit(comp, new_node1, new_node2):
    """
    Equivalente del coupler: un T per ciascuna linea + capacita' mutua tra i nodi centrali.
    nodes = [linea0 x=0, linea0 x=l, linea1 x=0, linea1 x=l] (stesso ordine di CPW_coupler.Y).
    """
    Cm = comp.C * comp.l                     # matrice di Maxwell per l'intera lunghezza
    C12 = abs(Cm[0, 1])
    C1g = max(Cm[0, 0] - C12, 1e-3 * Cm[0, 0])
    C2g = max(Cm[1, 1] - C12, 1e-3 * Cm[1, 1])
    L1 = comp.L[0, 0] * comp.l / 2
    L2 = comp.L[1, 1] * comp.l / 2
    n = comp.nodes
    return [L(n[0], new_node1, L1), C(new_node1, 0, C1g), L(new_node1, n[1], L1),
            L(n[2], new_node2, L2), C(new_node2, 0, C2g), L(new_node2, n[3], L2),
            C(new_node1, new_node2, C12)]


def _equivalent_netlist(netlist):
    """
    Netlist equivalente: toglie le resistenze e sostituisce CPW e coupler con i loro circuiti equivalenti.

    Returns
    -------
    equivalent_netlist : list
        Solo componenti a due terminali (C, L, J).
    mapping : dict
        id(componente originale) -> dict con gli indici dei rami equivalenti:
          C, L, J     : {'branch': i}
          CPW         : {'x_l': i_L1, 'x_0': i_L2, 'path': [i_L1, i_L2]}
          CPW_coupler : {'lines': [[i_L1, i_L2], [i_L3, i_L4]]}
    """
    reduced_netlist = _reduce_netlist(netlist) if _there_is_R(netlist) else netlist

    max_node = _max_node(reduced_netlist)
    equivalent_netlist = []
    mapping = {}
    for comp in reduced_netlist:
        start = len(equivalent_netlist)
        if isinstance(comp, CPW):
            max_node += 1
            equivalent_netlist.extend(_equivalent_CPW_circuit(comp, max_node))
            mapping[id(comp)] = {'x_l': start, 'x_0': start + 2, 'path': [start, start + 2]}
        elif isinstance(comp, CPW_coupler):
            max_node += 2
            equivalent_netlist.extend(_equivalent_CPW_coupler_circuit(comp, max_node - 1, max_node))
            mapping[id(comp)] = {'lines': [[start, start + 2], [start + 3, start + 5]]}
        else:
            equivalent_netlist.append(comp)
            mapping[id(comp)] = {'branch': start}
    return equivalent_netlist, mapping


# ---------------------------------------------------------------------------
# Matrici
# ---------------------------------------------------------------------------

def _node_index(eq_netlist):
    """Numerazione compatta dei nodi (il ground 0 resta 0). Evita colonne nulle dopo aver tolto le R."""
    nodes = sorted({n for c in eq_netlist for n in (c.node_minus, c.node_plus)} - {0})
    return {0: 0, **{n: i + 1 for i, n in enumerate(nodes)}}


def _Ced(eq_netlist):
    """
    Matrice delle capacita' nella base dei rami (diagonale).
    """
    return np.diag([c.C_value if isinstance(c, C) else 0.0 for c in eq_netlist])


def _W_matrix(eq_netlist, node_index):
    """
    Matrice di incidenza rami x nodi (senza ground): phi_ramo = phi_plus - phi_minus.
    """
    W = np.zeros((len(eq_netlist), len(node_index)))
    for i, c in enumerate(eq_netlist):
        W[i, node_index[c.node_minus]] -= 1
        W[i, node_index[c.node_plus]] += 1
    return W[:, 1:]


# ---------------------------------------------------------------------------
# Loop
# ---------------------------------------------------------------------------

def _oriented_loop(eq_netlist, branches, start_branch=None, start_node=None):
    """
    Vettore orientato (+1/-1/0 per ramo) di un loop semplice.

    Il loop viene percorso partendo da start_branch, uscendo da start_node.
    Un ramo percorso da node_minus a node_plus ha segno +1, altrimenti -1.
    Senza start_branch/start_node si parte dal primo ramo percorso da node_minus a node_plus
    (usato per i loop di completamento, che hanno flusso nullo e quindi verso irrilevante).
    """
    branches = list(dict.fromkeys(branches))
    deg = {}
    for b in branches:
        for n in (eq_netlist[b].node_minus, eq_netlist[b].node_plus):
            deg[n] = deg.get(n, 0) + 1
    if any(d != 2 for d in deg.values()):
        raise ValueError("The components of a loop must form a single closed path (each node touched exactly twice)")

    if start_branch is None:
        start_branch = branches[0]
        start_node = eq_netlist[start_branch].node_minus

    g = np.zeros(len(eq_netlist))
    c = eq_netlist[start_branch]
    if c.node_minus == start_node:
        g[start_branch], node = 1.0, c.node_plus
    else:
        g[start_branch], node = -1.0, c.node_minus
    end = start_node
    left = set(branches) - {start_branch}
    while left:
        nxt = next((b for b in left if node in (eq_netlist[b].node_minus, eq_netlist[b].node_plus)), None)
        if nxt is None:
            raise ValueError("The components of a loop do not form a single closed path")
        c = eq_netlist[nxt]
        if c.node_minus == node:
            g[nxt], node = 1.0, c.node_plus
        else:
            g[nxt], node = -1.0, c.node_minus
        left.remove(nxt)
    if node != end:
        raise ValueError("The loop is not closed")
    return g


def _loop_branches(loop_components, mapping):
    """
    Traduce i componenti di un loop (netlist originale) in rami della netlist equivalente
    e ricava il verso di percorrenza dall'ordine della lista.

    I componenti vanno elencati nell'ordine in cui si incontrano girando attorno al loop:
    il verso di percorrenza va dal primo componente al secondo. L'orientamento
    (node_minus/node_plus) dei singoli componenti non conta.
    Per un loop di due soli componenti (in parallelo) il verso e' quello del primo
    componente, da node_minus a node_plus.

    Una CPW e' percorsa lungo il conduttore centrale (i due induttori del T, non la capacita').
    Per un coupler si sceglie la linea i cui due estremi sono toccati dagli altri componenti del loop.

    Returns
    -------
    branches : list   rami del loop nella netlist equivalente
    start_branch : int  ramo da cui parte la percorrenza
    start_node : int    nodo da cui si esce percorrendo start_branch
    """
    nodes_others = set()
    for comp in loop_components:
        if hasattr(comp, "node_minus"):
            nodes_others |= {comp.node_minus, comp.node_plus}

    branches, comp_branches, comp_nodes = [], [], []
    for comp in loop_components:
        m = mapping.get(id(comp))
        if m is None:
            raise ValueError(f"Component {comp} of a loop is not in the (resistor-free) netlist")
        if 'branch' in m:
            br, nodes = [m['branch']], (comp.node_minus, comp.node_plus)
        elif 'path' in m:
            br, nodes = list(m['path']), (comp.node_minus, comp.node_plus)
        else:
            n = comp.nodes
            candidates = [k for k in range(2) if {n[2 * k], n[2 * k + 1]} <= nodes_others]
            if len(candidates) != 1:
                raise ValueError("Cannot determine which line of the CPW_coupler belongs to the loop")
            k = candidates[0]
            br, nodes = list(m['lines'][k]), (n[2 * k], n[2 * k + 1])
        branches.extend(br)
        comp_branches.append(br)
        comp_nodes.append(set(nodes))

    # verso di percorrenza dall'ordine della lista
    n_comp = len(loop_components)
    if n_comp >= 3:
        for i in range(n_comp):
            if not comp_nodes[i] & comp_nodes[(i + 1) % n_comp]:
                raise ValueError("List the components of a loop in the order you meet them going around the loop")
        shared = comp_nodes[0] & comp_nodes[1]
        if len(shared) != 1:
            raise ValueError("Cannot determine the loop direction from the first two components")
        exit_node = shared.pop()
        start_node = (comp_nodes[0] - {exit_node}).pop()
    else:
        first = loop_components[0]
        start_node = first.node_minus if hasattr(first, "node_minus") else min(comp_nodes[0])

    eq_first = comp_branches[0]
    # ramo del primo componente che tocca il nodo di partenza (per CPW/coupler: l'induttore a quell'estremo)
    return branches, eq_first, start_node


def _candidate_loops(eq_netlist):
    """
    Loop "corti" usati per completare la base: prima i loop di due rami in parallelo,
    poi la base di cicli minima del grafo semplice. Ognuno e' una lista di rami.
    """
    groups = {}
    for i, c in enumerate(eq_netlist):
        groups.setdefault(frozenset((c.node_minus, c.node_plus)), []).append(i)

    candidates = []
    G = nx.Graph()
    for key, idx in groups.items():
        if len(key) < 2:
            continue
        a, b = tuple(key)
        G.add_edge(a, b, branch=idx[0])
        for other in idx[1:]:
            candidates.append([idx[0], other])

    cycles = sorted(nx.minimum_cycle_basis(G), key=len)
    for cyc in cycles:
        sub = G.subgraph(cyc)
        if sub.number_of_edges() != len(cyc):
            continue
        candidates.append([d['branch'] for _, _, d in sub.edges(data=True)])
    for cyc in nx.cycle_basis(G):    # riserva, se la base minima non bastasse
        sub = G.subgraph(cyc)
        if sub.number_of_edges() == len(cyc):
            candidates.append([d['branch'] for _, _, d in sub.edges(data=True)])
    return candidates


def _as_pairs(loops):
    """
    Converte i loop in coppie (lista_di_componenti, flusso).
    Accetta oggetti con attributi .components e .flux (la classe loop di core.py)
    oppure tuple (componenti, flusso).
    """
    pairs = []
    for lp in loops:
        if hasattr(lp, "components") and hasattr(lp, "flux"):
            comps, flux = lp.components, lp.flux
        else:
            comps, flux = lp
        if len(comps) == 0:
            raise ValueError("A loop must contain at least one component")
        pairs.append((list(comps), float(flux)))
    return pairs


def _check_user_loops_are_meshes(eq_netlist, G_user, loops):
    """
    Avvisa se un loop dell'utente non e' una maglia (contiene altre maglie al suo interno).
    In quel caso il completamento automatico a flusso nullo non sa quale maglia interna
    contenga davvero il flusso, e la scelta diventa arbitraria.
    Il confronto e' fatto sulle coppie di nodi, quindi i rami in parallelo e la struttura
    interna delle CPW (sostituite dal loro T) non contano.
    """
    pairs_of = lambda branches: frozenset(frozenset((eq_netlist[b].node_minus, eq_netlist[b].node_plus))
                                          for b in branches)
    graph = nx.Graph()
    for c in eq_netlist:
        if c.node_minus != c.node_plus:
            graph.add_edge(c.node_minus, c.node_plus)
    meshes = set()
    for cyc in nx.minimum_cycle_basis(graph):
        sub = graph.subgraph(cyc)
        if sub.number_of_edges() == len(cyc):
            meshes.add(frozenset(frozenset(e) for e in sub.edges()))
    for row, (comps, _) in zip(G_user, loops):
        pairs = pairs_of(np.flatnonzero(row))
        if len(pairs) > 1 and pairs not in meshes:
            names = ", ".join(type(c).__name__ for c in comps)
            print(f"Warning: the loop [{names}] contains other meshes inside it. The flux is assigned to "
                  "one of them arbitrarily: define the loop as the mesh that actually contains the flux.")


def _G_matrix(eq_netlist, mapping, loops, n_loops):
    """
    Matrice dei vincoli di loop (n_loops x n_rami).

    Le prime righe sono i loop definiti dall'utente (nell'ordine dato); le altre completano
    la base con loop corti a flusso nullo (ipotesi: i loop dell'utente sono maglie
    fondamentali e tutte le altre maglie non concatenano flusso).
    """
    rows = []
    for comps, _ in loops:
        branches, first_branches, start_node = _loop_branches(comps, mapping)
        eq = eq_netlist
        start_branch = next(b for b in first_branches if start_node in (eq[b].node_minus, eq[b].node_plus))
        rows.append(_oriented_loop(eq_netlist, branches, start_branch, start_node))
    G = np.array(rows).reshape(len(rows), len(eq_netlist))
    _check_user_loops_are_meshes(eq_netlist, G, loops)
    if np.linalg.matrix_rank(G) < len(rows):
        raise ValueError("The user-defined loops are not independent")

    for cand in _candidate_loops(eq_netlist):
        if G.shape[0] == n_loops:
            break
        g = _oriented_loop(eq_netlist, cand)
        G_try = np.vstack([G, g])
        if np.linalg.matrix_rank(G_try) > G.shape[0]:
            G = G_try
    if G.shape[0] != n_loops:
        raise RuntimeError("Could not build a complete loop basis")

    # Con 3 o piu' rami in parallelo tra gli stessi due nodi, quali coppie di rami non
    # concatenano flusso dipende dal disegno: la scelta automatica e' arbitraria.
    groups = {}
    for i, c in enumerate(eq_netlist):
        groups.setdefault(frozenset((c.node_minus, c.node_plus)), []).append(i)
    big = [set(idx) for idx in groups.values() if len(idx) >= 3]
    for row in G[len(rows):]:
        nz = set(np.flatnonzero(row))
        if len(nz) == 2 and any(nz <= grp for grp in big):
            names = [f"{type(eq_netlist[k]).__name__}({eq_netlist[k].node_minus},{eq_netlist[k].node_plus})"
                     for k in sorted(nz)]
            print("Warning: ambiguous loop between parallel components " + ", ".join(names) +
                  " assumed to have zero flux. With 3 or more components in parallel, define the "
                  "zero-flux loops explicitly, e.g. loop([J, C_J], 0) and loop([L, C_L], 0).")
    return G


# ---------------------------------------------------------------------------
# Distribuzione irrotazionale
# ---------------------------------------------------------------------------

def irrotational_b(netlist, loops):
    """
    Coefficienti b del gauge irrotazionale.

    Parameters
    ----------
    netlist : list
        Netlist originale (con CPW, coupler, eventuali R).
    loops : list of loop (core.loop: .components, .flux) or of (list_of_components, phi_ext)
        Loop con flusso dipendente dal tempo. I componenti vanno elencati nell'ordine in cui
        si incontrano girando attorno al loop; un flusso positivo e' orientato con la regola
        della mano destra rispetto a quel verso (elenco in senso orario -> flusso entrante).

    Returns
    -------
    B : numpy.ndarray (n_rami_equivalenti x n_loops)
        B[k, e] = b del ramo k per un flusso unitario nel loop e.
    eq_netlist, mapping : la netlist equivalente e la mappa dai componenti originali.
    """
    loops = _as_pairs(loops)
    eq_netlist, mapping = _equivalent_netlist(netlist)
    node_index = _node_index(eq_netlist)
    W = _W_matrix(eq_netlist, node_index)
    nb, nn = W.shape

    graph = nx.MultiGraph()
    graph.add_nodes_from(node_index)
    graph.add_edges_from((c.node_minus, c.node_plus) for c in eq_netlist)
    n_loops = nb - len(node_index) + nx.number_connected_components(graph)

    # vincoli di nodo: sum_k C_k b_k = 0 su ogni nodo con almeno una capacita'
    K = W.T @ _Ced(eq_netlist)
    norms = np.abs(K).max(axis=1)
    if np.any(norms == 0):
        print("Warning: some nodes have no capacitance; the b distribution among their inductive branches "
              "is not unique (minimum-norm solution chosen, physically equivalent).")
    K = K[norms > 0] / norms[norms > 0, None]

    G = _G_matrix(eq_netlist, mapping, loops, n_loops)

    A = np.vstack([K, G])
    rhs = np.zeros((A.shape[0], len(loops)))
    rhs[K.shape[0]:K.shape[0] + len(loops), :] = np.eye(len(loops))
    B, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    if not np.allclose(A @ B, rhs, atol=1e-9):
        raise RuntimeError("The irrotational constraints have no solution")
    B[np.abs(B) < 1e-12] = 0.0
    return B, eq_netlist, mapping


def irrotational_phi_ext(netlist, loops):
    """
    Flussi esterni da assegnare ai componenti (gauge irrotazionale): phi_ext = sum_e b(e) * phi_ext(e).

    Returns
    -------
    dict : componente -> phi_ext
        J, L : numero
        CPW  : [fase in x=0 (node_plus), fase in x=l (node_minus)]
    I rami capacitivi non compaiono (il loro b non entra nell'Hamiltoniana).
    """
    loops = _as_pairs(loops)
    B, eq_netlist, mapping = irrotational_b(netlist, loops)
    phi = np.array([p for _, p in loops], dtype=float)
    b_tot = B @ phi

    result = {}
    for comp in netlist:
        m = mapping.get(id(comp))
        if m is None:
            continue
        if isinstance(comp, (J, L)):
            result[comp] = float(b_tot[m['branch']])
        elif isinstance(comp, CPW):
            result[comp] = [float(b_tot[m['x_0']]), float(b_tot[m['x_l']])]
        elif isinstance(comp, CPW_coupler):
            if np.any(np.abs(b_tot[sum(m['lines'], [])]) > 1e-12):
                raise NotImplementedError("Non-zero external flux on a CPW_coupler line is not supported")
    return result


def assign_irrotational_fluxes(netlist, loops):
    """
    Calcola i flussi esterni nel gauge irrotazionale e li assegna ai componenti della netlist
    di partenza (J, L, CPW), cosi' la simulazione dell'Hamiltoniana li usa direttamente.

    phi_ext finale = phi_ext statico assegnato dall'utente + sum_e b(e) * phi_ext(e).
    Il valore statico viene salvato in comp._phi_ext_static alla prima chiamata, quindi la
    funzione puo' essere richiamata (per esempio cambiando i valori dei flussi nei loop)
    senza sommare due volte il contributo dinamico.

    Parameters
    ----------
    netlist : list
        La netlist di partenza (la stessa passata a QCircuit: i componenti sono gli stessi oggetti).
    loops : list of loop (core.loop: .components, .flux) or of (list_of_components, phi_ext)
        Loop con flusso dipendente dal tempo.
    """
    dynamic = irrotational_phi_ext(netlist, loops)

    for comp in netlist:
        if not isinstance(comp, (J, L, CPW)):
            continue
        if not hasattr(comp, "_phi_ext_static"):
            comp._phi_ext_static = deepcopy(comp.phi_ext)
        static = comp._phi_ext_static
        dyn = dynamic.get(comp)

        if isinstance(comp, CPW):
            static_arr = np.zeros(2) if len(static) == 0 else np.asarray(static, dtype=float)
            dyn_arr = np.zeros(2) if dyn is None else np.asarray(dyn, dtype=float)
            total = static_arr + dyn_arr
            comp.phi_ext = [] if (len(static) == 0 and not np.any(total != 0)) else [float(x) for x in total]
        else:
            comp.phi_ext = float(static) + (0.0 if dyn is None else dyn)