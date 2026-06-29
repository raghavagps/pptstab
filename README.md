# **PPTStab**
Designing of thermostable proteins with a desired melting temperature

## Introduction
PPTStab is developed to predict the thermostability of proteins and design thermostable proteins. The standalone version uses an ANN+MLP ensemble regressor model.
PPTStab is also available as a web-server at https://webs.iiitd.edu.in/raghava/pptstab. Please read/cite the content about PPTStab for complete information including the algorithm behind the approach.

## Installation
To use PPTStab, follow these steps to set up your environment:

### Install using environment.yml

1. Create a new Conda environment from the environment.yml file, replacing `<env_name>` with your preferred name (e.g. `pptstab`):

```bash
conda env create -n <env_name> -f environment.yml
```

Example:

```bash
conda env create -n pptstab -f environment.yml
```

> **Note:** If an environment with that name already exists, remove it first:
> ```bash
> conda env remove -n <env_name>
> conda env create -n <env_name> -f environment.yml
> ```

2. Activate the environment:

```bash
conda activate pptstab
```

### Required Dependencies
The standalone version of PPTStab is written in Python 3.11. The following libraries are required:

| Package | Version |
|---|---|
| numpy | 1.26.4 |
| scikit-learn | 1.6.1 |
| tensorflow | 2.16.1 |
| torch | 2.3.1 |
| transformers | 4.41.2 |
| pandas | 2.3.3 |
| onnxruntime | 1.18.0 |
| joblib | 1.4.2 |
| tqdm | 4.66.4 |

## Minimum Usage
To see all available options, type:

```bash
python pptstab.py -h
```

To run the example (predict melting temperature using default settings):

```bash
python pptstab.py -i example.fasta -f 1
```

Results will be saved to `outfile.csv` by default.

## Full Usage

```
usage: pptstab.py [-h] [-i INPUT] [-o OUTPUT] [-j {1,2}] [-d {1,2}] [-f {0,1}] [-m {EMB,AAC,SER}]
```

### Arguments

| Argument | Description | Default |
|---|---|---|
| `-i INPUT` | Input protein/peptide sequence(s) in FASTA format or one sequence per line | Required |
| `-o OUTPUT` | Output file name | outfile.csv |
| `-j {1,2}` | Job type: 1 = Predict, 2 = Design | 1 |
| `-f {0,1}` | Flag: 1 = lysate, 0 = cell | 1 |
| `-d {1,2}` | Display: 1 = thermophilic proteins only, 2 = all proteins | 1 |
| `-m {EMB,AAC,SER}` | Method: EMB (ProtBERT embeddings), AAC (amino acid composition), SER (Shannon entropy) | EMB |

### Usage Examples

**Predict Tm using default ProtBERT embeddings:**
```bash
python pptstab.py -i example.fasta -f 1 -m EMB
```

**Predict Tm using SER method (lysate context):**
```bash
python pptstab.py -i example.fasta -f 1
```

**Predict Tm using AAC method, cell context, show all results:**
```bash
python pptstab.py -i example.fasta -f 0 -m AAC -d 2
```

**Design thermostable mutants using SER method:**
```bash
python pptstab.py -i example.fasta -j 2 -f 1 -m SER -o design_output.csv
```

**Input File:** Accepts sequences in FASTA format or one sequence per line (single-letter amino acid code).

**Output File:** Results saved in CSV format. If no output file is specified, results are stored in `outfile.csv`.

**Job:** 1 = Predict melting temperature (Tm) for input sequences. 2 = Design mode — generates all single-point mutants of input sequences and predicts Tm for each.

**Flag:** Sets the experimental context. `1` = lysate (cell-free system), `0` = cell. Default is `1` (lysate).

**Method:** Selects the feature/model type. `SER` uses Shannon entropy per residue, `AAC` uses amino acid composition, `EMB` uses ProtBERT protein embeddings (requires GPU for large inputs).

## PPTStab Package Files

| File | Description |
|---|---|
| `INSTALLATION` | Installation instructions |
| `LICENSE` | License information |
| `README.md` | This file |
| `pptstab.py` | Main Python program |
| `environment.yml` | Conda environment file |
| `example.fasta` | Example input file with protein sequences in FASTA format |
| `example_predict_output.csv` | Example output for predict module |
| `example_design_output.csv` | Example output for design module |
| `models/` | Pre-trained model files (ANN ONNX + MLP pkl) for each method |
