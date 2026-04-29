# **PPTStab**
Designing of thermostable proteins with a desired melting temperature
## Introduction
PPTStab is developed to predict the thermostability of proteins and design the thermostable proteins. In the standalone version, ANN+MLP ensemble regressor based model.
PPTStab is also available as web-server at https://webs.iiitd.edu.in/raghava/pptstab. Please read/cite the content about the PPTStab for complete information including algorithm behind the approach.

## Publication
Tijare, P., Kumar, N. & Raghava, G.P.S. Prediction and design of thermostable proteins with a desired melting temperature. Sci Rep 15, 16683 (2025). https://doi.org/10.1038/s41598-025-98667-9


## Installation
To use pptstab, follow these steps to set up your environment:

## Install using environment.yml

1. Create a new Conda environment from the environment.yml file:

```
conda env create -f environment.yml
```
2. Activate the newly created environment:
```
conda activate pptstab
```

## Standalone
The Standalone version of PPTStab is written in python3 and following libraries are necessary for the successful run:
- scikit-learn==1.0.2
- transformers==4.44.2
- tensorflow==2.13.0
- Pandas==2.0.3
- Numpy==1.22.4
- torch==2.4.1

## Minimum USAGE
To know about the available option for the stanadlone, type the following command:
```
python pptstab.py -h
```
To run the example, type the following command:
```
python3 pptstab.py -i example_input.fa -f 1
```
This will predict if the submitted sequences melting temperature (Tm). It will use other parameters by default. It will save the output in "outfile.csv" in CSV (comma seperated variables).

## Full Usage
```
usage: pptstab.py [-h] 
                  [-i INPUT]
                  [-o OUTPUT]
                 [-j {1,2}]
                 [-d {1,2}]
    	         [-f {0,1}]
                 [-m {EMB,AAC,SER}]

```
```
Please provide following arguments for successful run

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence(s) in FASTA format or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -j {1,2}, --job {1,2}
                        Job Type: 1:Predict, 2: Design, by default 1
  -f FLAG, --flag FLAG {0,1} 
                        Cell Flag: Value between 0 or 1 by default 1
  -d {1,2}, --display {1,2}
                        Display: 1:Thermophilic proteins only, 2: All peptides, by default 1
  -m {EMB,AAC,SER}, --method {EMB,AAC,SER}
                        Display: EMB for using the embedding model (ProtBert), SER for Shannon entropy for all residues model, AAC for amino acid composition model, by default SER

```

**Input File:** It allow users to provide input in the FASTA format.

**Output File:** Program will save the results in the CSV format, in case user do not provide output file name, it will be stored in "outfile.csv".

**Job:** User is allowed to choose between three different modules, such as, 1 for prediction, and 2 for Designing, by default its 1.

**flag**: User can set the ‘lysate’ or ‘cell’ flag.

**method:** This option allow users to select the model based on composition and embeddings, by default SER.

PPTStab Package Files
=======================
It contantain following files, brief descript of these files given below

INSTALLATION                    : Installations instructions

LICENSE                         : License information

README.md                       : This file provide information about this package

pptstab.py                      : Main python program

example.fasta                   : Example file contain peptide sequenaces in FASTA format

example_predict_output.csv      : Example output file for predict module

example_design_output.csv       : Example output file for design module
