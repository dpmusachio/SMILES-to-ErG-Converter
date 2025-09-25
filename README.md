# SMILES → ErG Converter

A Python pipeline to convert molecules from SMILES into **Extended Reduced Graphs (ErGs)**.  
Nodes represent pharmacophoric features (donors, acceptors, charges, hydrophobes, aromatics),  
edges encode shortest-path distances.  

Based on:  
M. Stiefl, M. Baumann, K. Hähnke, G. Zaliani.  
*ErG: 2D Pharmacophore Descriptions for Scaffold Hopping.*  
**Journal of Chemical Information and Modeling** 2006, 46 (1), 208–220.  
[https://doi.org/10.1021/ci050344m](https://doi.org/10.1021/ci050344m)

## ⚙️ Installation

Create a conda environment with the provided `environment.yml` file and activate it.


conda env create -f environment.yml


conda activate erg



Dependencies:
- Python 3.10  
- RDKit  
- NetworkX  

Once in the conda environment, run:

python src/smiles_to_erg.py input/ligands.smi output/ligands.erg.jsonl

with whichever input you want

## 📖 Reference

M. Stiefl, M. Baumann, K. Hähnke, G. Zaliani.  
*ErG: 2D Pharmacophore Descriptions for Scaffold Hopping.*  
**J. Chem. Inf. Model.** 2006, 46 (1), 208–220.  
[https://doi.org/10.1021/ci050344m](https://doi.org/10.1021/ci050344m)

---

## 📜 License

MIT License.


