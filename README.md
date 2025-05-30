# 🧬 DeepPurpose-MD-Discovery

A streamlined, Colab-optimized drug discovery pipeline integrating:

- ✅ Ligand-target prediction with a custom-trained [DeepPurpose](https://github.com/mosmos6/Deeppurpose) fork
- ✅ Structural docking using [AutoDock Vina](http://vina.scripps.edu/)
- ✅ GPU-accelerated Molecular Dynamics with [OpenMM](https://openmm.org/) and [OpenFF](https://openforcefield.org/)
- ✅ Stability and mechanistic analyses via PCA, FEL, RMSD, H-bonding, water networks, and more.

This repo demonstrates how to simulate and evaluate ligand–protein interactions from end-to-end using SARS-CoV-2 viral proteins as case studies.

---

## 🔧 Setup Instructions (Google Colab)

This pipeline is designed for use in **Google Colab**, with full support for `condacolab`. A demo notebook is included in this repository to reproduce all steps.

---

### ✅ Step 1: Enable Conda in Colab

Paste the following **at the very top of your Colab notebook**:

```python
!pip install -q condacolab
import condacolab
condacolab.install()
```



🔄 NOTE: This will crash your runtime once. That's expected.

After Colab restarts, rerun the following cell:

```python
import condacolab
condacolab.check()
```

### ✅ Step 2: Clone Required Repositories

```python
# Main pipeline repo (this one)
!git clone https://github.com/BioMolDynamics/DeepPurpose-MD-Discovery.git

# Custom fork of DeepPurpose (installed later)
!git clone https://github.com/BioMolDynamics/Deeppurpose
```

### ✅ Step 3: Download AutoDock Vina Binary

```python
!wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
!chmod +x vina_1.2.5_linux_x86_64
!./vina_1.2.5_linux_x86_64 --version
```

### ✅ Step 4: Install Conda Environment

```python
!mamba env create -f environment.yml
```

### ✅ Step 5: Finalize Setup

```python
# Install custom fork of DeepPurpose without overwriting key dependencies
!conda run -n deeppurpose-md-env pip install --no-deps ./Deeppurpose

# Install optional dependencies like ProtTrans
!conda run -n deeppurpose-md-env python scripts/install_optional.py
```

### ✅ Step 6: Run the Pipeline

```python
!conda run -n deeppurpose-md-env python scripts/1_prepare_ligand.py "$ligand_smiles"
!conda run -n deeppurpose-md-env python scripts/2_prepare_receptor.py "$pdb_id" --strict-protein
!conda run -n deeppurpose-md-env python scripts/3_docking_vina.py --use-residue-centroid
!conda run -n deeppurpose-md-env python scripts/3b_prepare_protein.py
!conda run -n deeppurpose-md-env python scripts/4_align_ligand.py
!conda run -n deeppurpose-md-env python scripts/5_md_simulation.py
!conda run -n deeppurpose-md-env python scripts/5b_md_analysis_RMSD_RMSF.py
!conda run -n deeppurpose-md-env python scripts/6_md_analysis.py
!conda run -n deeppurpose-md-env python scripts/7_deeppurpose_training.py
!conda run -n deeppurpose-md-env python scripts/8_deeppurpose_prediction.py
```
Each script corresponds to a specific stage in the full drug discovery pipeline — from ligand design to MD simulation to deep learning prediction.


### 🧪 Dataset Information

This demo uses a COVID-19 specific subset of BindingDB, available from UC San Diego.  
To simplify setup, we provide pre-cleaned versions of this dataset:

- `BindingDB_Covid-19.tsv` (214MB, hosted via OSF)
- `strong_binders_cleaned.csv` (optional for filtering)
- `protein.faa` (optional for filtering)
- `metrics - SARS2 FASTA.csv` (matching data of SARS-CoV-2 proteins and FASTA)

📎 Dataset download link: [https://osf.io/your-dataset-id](https://osf.io/your-dataset-id)

You are welcome to use your own SMILES/FASTA data by modifying `7_deeppurpose_training.py`.

### 🧬 Features
Compatible with custom ligands and proteins (via CLI args or code edit)

Supports residue-based centroid docking

Robust fallback handling for PDBFixer and Meeko failures

Cleaned outputs ready for visualization and MD refinement

### 📜 License
MIT License. Please cite this repository if used in academic work.

### 📫 Contact
Maintained by BioMolDynamics
For academic inquiries, collaboration, or feedback, please open an issue.
