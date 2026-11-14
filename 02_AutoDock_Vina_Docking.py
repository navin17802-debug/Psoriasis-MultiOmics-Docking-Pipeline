# ==============================================================================
# MASTER SCRIPT 2: AUTOMATED AUTODOCK VINA MOLECULAR DOCKING (PYTHON)
# Evaluates 5 Hub Genes x 8 Botanical Oil Phytoconstituents (40 Complexes)
# ==============================================================================

import os
import re
import shutil
import subprocess
import urllib.request
import numpy as np
import pandas as pd

# 1. Install Vina and OpenBabel if running on Linux / Google Colab
print(">>> [1/5] Checking & Installing AutoDock Vina and OpenBabel...")
subprocess.run(
    "apt-get update -qq && apt-get install -y -qq autodock-vina openbabel",
    shell=True,
    check=True
)

# 2. Setup clean directory workspace
folders = [
    "docking_workspace/proteins_raw",
    "docking_workspace/proteins_pdbqt",
    "docking_workspace/ligands_sdf",
    "docking_workspace/ligands_pdbqt",
    "docking_workspace/docked_complexes",
    "docking_workspace/final_results"
]
for folder in folders:
    os.makedirs(folder, exist_ok=True)

# 3. Define Target Proteins (5 Hub Genes) and Ligands (8 Phytoconstituents)
targets = {
    "IL1B": "5R8Q",
    "STAT1": "1BF5",
    "TLR2": "2Z7X",
    "CXCL8": "5D14",
    "CXCL10": "1O7Y"
}

ligands = {
    "Alpha_Linolenic_Acid": {"cid": "5280934", "smiles": None},
    "Gamma_Tocopherol":     {"cid": "9229",    "smiles": None},
    "Punicic_Acid":         {"cid": "5281126", "smiles": None},
    "Beta_Sitosterol":      {"cid": "222284",  "smiles": None},
    "Alpha_Pinene":         {"cid": "6654",    "smiles": None},
    "Cedrol":               {"cid": "65575",   "smiles": None},
    "Oleic_Acid":           {"cid": "445639",  "smiles": None},
    "Niazimicin":           {
        "cid": "5471459",
        "smiles": "CCOC(=S)NCC1=CC=C(C=C1)OC2C(C(C(C(O2)C)O)O)O"
    }
}

STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}

# 4. Download, Clean, and Convert 5 Target Proteins to PDBQT
print("\n>>> [2/5] Preparing 5 Hub Gene Receptor Structures (PDB -> PDBQT)...")
grid_centers = {}

for gene, pdb_id in targets.items():
    raw_pdb = f"docking_workspace/proteins_raw/{gene}_{pdb_id}_raw.pdb"
    clean_pdb = f"docking_workspace/proteins_raw/{gene}_{pdb_id}_clean.pdb"
    receptor_pdbqt = f"docking_workspace/proteins_pdbqt/{gene}.pdbqt"

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, raw_pdb)

    selected_chain = None
    coords = []
    with open(raw_pdb, "r") as fin, open(clean_pdb, "w") as fout:
        for line in fin:
            if line.startswith("ATOM"):
                alt_loc = line[16]
                res_name = line[17:20].strip()
                chain_id = line[21]
                if res_name not in STD_AA or alt_loc not in (" ", "A"):
                    continue
                if selected_chain is None:
                    selected_chain = chain_id
                if chain_id == selected_chain:
                    fout.write(line)
                    coords.append([
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54])
                    ])
        fout.write("END\n")

    center = np.array(coords).mean(axis=0)
    grid_centers[gene] = (round(center[0], 3), round(center[1], 3), round(center[2], 3))

    cmd = f"obabel {clean_pdb} -O {receptor_pdbqt} -xr -p 7.4 --partialcharge gasteiger"
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    print(f"    [OK] {gene} ({pdb_id}) | Grid Center: {grid_centers[gene]}")

# 5. Download, Minimize (MMFF94), and Convert 8 Ligands to PDBQT (Zero-NaN Guaranteed)
print("\n>>> [3/5] Preparing 8 Phytoconstituent Ligands (MMFF94 3D -> PDBQT)...")
for lig_name, info in ligands.items():
    cid = info["cid"]
    smiles = info["smiles"]
    sdf_file = f"docking_workspace/ligands_sdf/{lig_name}.sdf"
    lig_pdbqt = f"docking_workspace/ligands_pdbqt/{lig_name}.pdbqt"

    if smiles is not None:
        # Direct 3D build from SMILES for Niazimicin to avoid PubChem 2D NaN issue
        cmd_3d = f'obabel -:"{smiles}" -O {sdf_file} --gen3d -h --minimize --ff MMFF94 --steps 500'
        subprocess.run(cmd_3d, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        cmd_pdbqt = f"obabel {sdf_file} -O {lig_pdbqt} --partialcharge gasteiger"
        subprocess.run(cmd_pdbqt, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    else:
        url_3d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        urllib.request.urlretrieve(url_3d, sdf_file)
        cmd_pdbqt = f"obabel {sdf_file} -O {lig_pdbqt} --minimize --ff MMFF94 --steps 250 --partialcharge gasteiger"
        subprocess.run(cmd_pdbqt, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    print(f"    [OK] {lig_name} (CID: {cid}) converted to 3D PDBQT")

# 6. Run AutoDock Vina Across All 40 Combinations
print("\n>>> [4/5] Running AutoDock Vina Simulations (40 Protein-Ligand Pairs)...")
results_list = []
matrix_data = {gene: {} for gene in targets.keys()}
current_run = 0
total_runs = len(targets) * len(ligands)

for gene, pdb_id in targets.items():
    receptor_pdbqt = f"docking_workspace/proteins_pdbqt/{gene}.pdbqt"
    cx, cy, cz = grid_centers[gene]

    for lig_name, info in ligands.items():
        current_run += 1
        cid = info["cid"]
        lig_pdbqt = f"docking_workspace/ligands_pdbqt/{lig_name}.pdbqt"
        out_pose = f"docking_workspace/docked_complexes/{gene}_{lig_name}_docked.pdbqt"

        vina_cmd = [
            "vina",
            "--receptor", receptor_pdbqt,
            "--ligand", lig_pdbqt,
            "--center_x", str(cx),
            "--center_y", str(cy),
            "--center_z", str(cz),
            "--size_x", "30",
            "--size_y", "30",
            "--size_z", "30",
            "--exhaustiveness", "8",
            "--out", out_pose
        ]

        proc = subprocess.run(vina_cmd, capture_output=True, text=True, check=True)

        best_score = np.nan
        if os.path.exists(out_pose):
            with open(out_pose, "r") as f_out:
                for line in f_out:
                    if "REMARK VINA RESULT:" in line:
                        best_score = float(line.split()[3])
                        break
        if np.isnan(best_score):
            match = re.search(r"^\s*1\s+(-?\d+\.\d+)", proc.stdout, re.MULTILINE)
            if match:
                best_score = float(match.group(1))

        matrix_data[gene][lig_name] = best_score
        results_list.append({
            "Hub_Gene": gene,
            "PDB_ID": pdb_id,
            "Phytoconstituent_Ligand": lig_name.replace("_", " "),
            "PubChem_CID": cid,
            "Binding_Affinity_kcal_mol": best_score
        })
        print(f"    [{current_run:02d}/{total_runs}] {gene:6s} + {lig_name:20s} --> {best_score:.3f} kcal/mol")

# 7. Export Clean CSV Tables for R Visualization
print("\n>>> [5/5] Exporting Final CSV Tables for R Pipeline...")
df_long = pd.DataFrame(results_list)
df_long.to_csv("docking_workspace/final_results/Full_40_Docking_Scores_Long.csv", index=False)

df_matrix = pd.DataFrame(matrix_data).T
df_matrix.columns = [c.replace("_", " ") for c in df_matrix.columns]
df_matrix["Mean_Affinity (kcal/mol)"] = df_matrix.mean(axis=1).round(2)
df_matrix["Best_Affinity (kcal/mol)"] = df_matrix.min(axis=1).round(2)
df_matrix = df_matrix.sort_values(by="Mean_Affinity (kcal/mol)", ascending=True)
df_matrix.to_csv("docking_workspace/final_results/Docking_Affinity_Matrix_Ranked.csv")

# Also save a copy in root folder so R script finds it immediately
df_matrix.to_csv("Docking_Affinity_Matrix_Ranked.csv")

print("\n✅ Python Docking Pipeline Complete! CSV saved as 'Docking_Affinity_Matrix_Ranked.csv'.")
