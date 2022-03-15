code of the paper 'SCP4SSD: a Serverless Cloud Platform for the prediction of nucleotide Sequence Synthesis Difficulty'

## Requirements

- Python=3.9
- biopython=1.78
- pandas=1.2.4
- scikit-learn=0.24.2
- numpy=1.20.2
- auto-sklearn

## Installation

1. Clone the repo

```shell
git clone https://github.com/JustinDoIt/scp4ssd.git
cd scp4ssd
```

2. Create Anaconda Environment

```shell
conda env create -f name.yml
```

3. Activate the environment

```shell
conda activate scp4ssd
```

## Usage & Example

```shell
python predict.py --fasta example.fna --out example_out.csv
```

## Cite us

If this repo help you, happy to cite our paper (coming soon...)

## License

Distributed under the MIT License. See LICENSE for more information
