import os
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig

fdef = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
feat_factory = ChemicalFeatures.BuildFeatureFactory(fdef)

def assign_features(mol):
    feats = feat_factory.GetFeaturesForMol(mol)
    donors, acceptors, pos, neg = set(), set(), set(), set()
    for f in feats:
        fam = f.GetFamily().lower()
        if fam == 'donor': donors.update(f.GetAtomIds())
        elif fam == 'acceptor': acceptors.update(f.GetAtomIds())
    for a in mol.GetAtoms():
        ch = a.GetFormalCharge()
        if ch > 0: pos.add(a.GetIdx())
        if ch < 0: neg.add(a.GetIdx())
    flipflop = donors & acceptors
    return donors, acceptors, pos, neg, flipflop

def find_endcaps(mol):
    endcaps, thio = set(), set()
    # TODO: refine SMARTS; placeholder implementation
    return endcaps, thio

def ring_centroids_and_keep_atoms(mol):
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    systems = [set(r) for r in atom_rings]
    centroids = []
    keep_atoms = set()
    ring_atoms = {a for r in systems for a in r}

    for i, sys_atoms in enumerate(systems):
        aromatic_flag = any(mol.GetAtomWithIdx(j).GetIsAromatic() for j in sys_atoms)
        centroids.append({'type':'Ar' if aromatic_flag else 'Hf', 'sys_atoms':list(sys_atoms)})
        # keep substituted atoms
        for j in sys_atoms:
            a = mol.GetAtomWithIdx(j)
            if any(n.GetIdx() not in sys_atoms for n in a.GetNeighbors()):
                keep_atoms.add(j)

    return centroids, keep_atoms, ring_atoms
