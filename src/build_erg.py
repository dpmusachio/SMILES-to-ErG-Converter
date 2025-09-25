import networkx as nx
from utils import assign_features, find_endcaps, ring_centroids_and_keep_atoms

def build_erg_graph(mol):
    donors, acceptors, pos, neg, flipflop = assign_features(mol)
    endcaps, thio = find_endcaps(mol)
    centroids, keep_ring_atoms, ring_atoms = ring_centroids_and_keep_atoms(mol)

    # Base graph for distances
    Gbase = nx.Graph()
    kept_atoms = {a.GetIdx() for a in mol.GetAtoms() if not a.IsInRing()} | keep_ring_atoms
    for a in kept_atoms:
        Gbase.add_node(a)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in kept_atoms and j in kept_atoms:
            Gbase.add_edge(i, j)

    centroid_map = {}
    for k, c in enumerate(centroids):
        cid = f"CENT{k}"
        centroid_map[cid] = c['sys_atoms']
        Gbase.add_node(cid)
        for i in c['sys_atoms']:
            if i in keep_ring_atoms:
                Gbase.add_edge(cid, i)

    PPs, rep_map = [], {}
    # ring centroids
    for cid, atoms in centroid_map.items():
        PPs.append((cid, c['type'], {'origin': 'ring_centroid', 'atoms': atoms}))
        rep_map[cid] = cid

    # donors/acceptors/charges
    for i in range(mol.GetNumAtoms()):
        if i in donors:
            PPs.append((f"D_{i}", 'D', {'origin':'atom','atoms':[i]}))
            rep_map[f"D_{i}"] = i
        if i in acceptors:
            PPs.append((f"A_{i}", 'A', {'origin':'atom','atoms':[i]}))
            rep_map[f"A_{i}"] = i
        if i in pos:
            PPs.append((f"P_{i}", '+', {'origin':'atom','atoms':[i]}))
            rep_map[f"P_{i}"] = i
        if i in neg:
            PPs.append((f"N_{i}", '-', {'origin':'atom','atoms':[i]}))
            rep_map[f"N_{i}"] = i
    for i in sorted(endcaps | thio):
        PPs.append((f"E_{i}", 'Hf', {'origin':'endcap','atoms':[i]}))
        rep_map[f"E_{i}"] = i

    erg = nx.Graph()
    for nid, ntype, data in PPs:
        erg.add_node(nid, type=ntype, **data)

    nlist = [nid for nid,_,_ in PPs]
    for u_idx in range(len(nlist)):
        for v_idx in range(u_idx+1, len(nlist)):
            u, v = nlist[u_idx], nlist[v_idx]
            try:
                rep_u, rep_v = rep_map[u], rep_map[v]
                d = nx.shortest_path_length(Gbase, rep_u, rep_v)
                erg.add_edge(u, v, distance=int(d))
            except Exception:
                pass

    return erg
