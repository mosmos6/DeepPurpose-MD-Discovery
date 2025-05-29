# 🧬 DeepPurpose-MD-Discovery

A streamlined, Colab-optimized drug discovery pipeline integrating:

- ✅ Ligand-target prediction with a custom-trained [DeepPurpose](https://github.com/mosmos6/Deeppurpose) fork
- ✅ Structural docking using [AutoDock Vina](http://vina.scripps.edu/)
- ✅ GPU-accelerated Molecular Dynamics with [OpenMM](https://openmm.org/) and [OpenFF](https://openforcefield.org/)
- ✅ Stability and mechanistic analyses via PCA, FEL, RMSD, H-bonding, water networks, and more.

This repo demonstrates how to simulate and evaluate ligand–protein interactions from end-to-end using SARS-CoV-2 viral proteins as case studies.

---

## 🔧 Setup Instructions (Google Colab)

This pipeline is designed for use in **Google Colab**, with full support for `condacolab`.

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
!mamba env update -n base -f environment.yml
```

### ✅ Step 2: Download AutoDock Vina Binary

AutoDock Vina is not available via pip/conda, so download it manually:

```python
!wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
!chmod +x vina_1.2.5_linux_x86_64
!./vina_1.2.5_linux_x86_64 --version
```

### ✅ Step 3: Clone & Install Custom DeepPurpose Fork

```python
!git clone https://github.com/BioMolDynamics/Deeppurpose
!pip install ./Deeppurpose
```

### ✅ Step 4: Run the Pipeline

You can now run the full discovery pipeline:

- DeepPurpose prediction with ProtTrans embeddings

- Ligand docking with Vina

- MD simulation using OpenMM/OpenFF

- Output analysis (PCA, FEL, H-bond, RMSD, etc.)

Sample .py and .ipynb scripts are included in the repository.

### 📂 Repository Contents

environment.yml — Full environment for conda setup

pipeline_main.py — Unified script for ligand → DeepPurpose → docking → MD → analysis

notebooks/ — Example Colab notebooks (coming soon)

docs/ — Description of each pipeline module (WIP)

### 📜 License
MIT License. Please cite this repository if used in academic work.
