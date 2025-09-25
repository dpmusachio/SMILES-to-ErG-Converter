import sys, json
from rdkit import Chem
from build_erg import build_erg_graph

def smiles_to_erg_jsonl(smiles_file, out_file):
    with open(smiles_file, "r") as fin, open(out_file, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            parts = line.strip().split()
            smi = parts[0]
            mol_id = parts[1] if len(parts) > 1 else None

            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue

            erg = build_erg_graph(mol)

            nodes = [
                {"id": n, "type": d.get("type"),
                 "atoms": d.get("atoms", []),
                 "origin": d.get("origin")}
                for n, d in erg.nodes(data=True)
            ]
            edges = [
                {"u": u, "v": v, "distance": d["distance"]}
                for u, v, d in erg.edges(data=True)
            ]

            fout.write(json.dumps({
                "id": mol_id,
                "smiles": smi,
                "nodes": nodes,
                "edges": edges
            }) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python src/smiles_to_erg.py input.smi output.jsonl")
        sys.exit(1)
    smiles_to_erg_jsonl(sys.argv[1], sys.argv[2])
