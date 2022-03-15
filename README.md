code of the paper 'SCP4SSD: a Serverless Cloud Platform for the prediction of nucleotide Sequence Synthesis Difficulty'

## Requirements

- Python=3.9
- biopython=1.78
- pandas=1.2.4
- scikit-learn=0.24.2
- numpy=1.20.2

## Install

1. clone the repo

```shell
git clone https://github.com/JustinDoIt/scp4ssd.git
```

2. create the conda env

```shell
conda env create -f name.yml
```

or install python packages via pip

```shell
cd scp4ssd
pip install -r requirements.txt
```

## Usage & Example

```shell
cd scp4ssd
python predict.py --fasta example.fna --out example_out.csv
```

## Cite us

If this repo help you, happy to cite our paper (coming soon...)

## License

Distributed under the MIT License. See LICENSE for more information
