# Nanopore sequencing FP-chimera-classifier
This repository provides a deep learning-based filtering method for identifying and removing false positive (FP) chimeras generated during Nanopore sequencing, specifically for data prepared using the ligation method. The tool leverages a ResNet-based model to improve the accuracy of downstream genomic analyses by effectively filtering out chimeric reads.

## Installation
To set up the FP-chimera-classifier, follow the steps below. We recommend using ```conda``` to manage the computing environment.
```bash
# Clone the Repository
git clone https://github.com/xzhbio/FP-chimera-classifier.git

# Navigate to the Directory
cd FP-SVs-classifier

# Set Up the Environment
conda env create -n FP-chimera-classifier_env -f env.yaml
```

## Usage

The `main.py` script is the core of the FP-chimera-classifier. Below are the required and optional arguments for running the tool.

### Command-Line Options

| Argument          | Description                                                                 | Required/Optional |
|-------------------|-----------------------------------------------------------------------------|-------------------|
| `--input`, `-i`   | Path to the input FASTQ file.                                               | **Required**      |
| `--ref`           | Path to the reference genome file.                                          | **Required**      |
| `--fast5`         | Path to the directory containing FAST5 files.                               | **Required**      |
| `--output`, `-o`  | Path to the output directory (default: `./`).                               | Optional          |
| `--cuda`          | Specify the GPU(s) to use (e.g., `0` for GPU 0, `0,1,2` for multiple GPUs). | Optional          |
| `--chimeric`      | If provided, outputs potential cross- and within-chimera PAF files.         | Optional          |
| `--model`         | Choose the pre-trained model: `NA12878` or `microbe`.                       | **Required**      |

### Example Command

```bash
python main.py --input sample.fastq --ref reference.fasta --fast5 fast5_directory --output results --cuda 0 --model NA12878
```

###  Pre-Trained Models
The FP-chimera-classifier provides two pre-trained models:
  1.NA12878: Optimized for human genomic data.
  2.Microbe: Designed for microbial genomic data.
Choose the appropriate model based on your dataset.

## Ootput
The tool generates the following outputs:

- **Filtered Reads**: A txt file containing FP chimeric reads id.
- **Chimera PAF Files** (optional): If the `--chimeric` flag is provided, the tool outputs PAF files for potential cross- and within-chimeras.
