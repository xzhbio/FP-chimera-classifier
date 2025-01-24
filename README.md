# Nanopore sequencing FP-chimera-classifier
In this repository, we provide a filter method for the FP chimera generated through nanopore sequencing which is using the ligation preparation method.

# Installation
In this part, we introcduce the installation step for FP-chimera-classifier. We use ```conda``` to manage the computing environment.
```
# download the source code
git clone https://github.com/xzhbio/FP-chimera-classifier.git

# enter the directory
cd FP-SVs-classifier

# install dependencies. Only test on Ubuntu System
conda env create -n FP-chimera-classifier_env -f env.yaml
```

# Usage
The main.py requires several options:
  1.required: the path of FASTQ file.
  2.required: the path of reference file.
  3.required: the path of FAST5 directory
  4.optional: the path of output directory, default is ```./```
  5.optional: choose the GPU for process, default is ```0``` (e.g. ```0,1,2``` for 3 GPU)
  6.optional: choose whether to output potential FP cross and within chimera PAF file.
  7.required: choose the model used for classifier, ```NA12878``` or ```microbe```.
```
usage: main.py [-h] --input INPUT --ref REF --fast5 FAST5 [--output OUTPUT] [--cuda CUDA] [--chimeric] [--model {NA12878,microbe}]

FP-SVs classifier

options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        input path of fastq file
  --ref REF             reference path
  --fast5 FAST5         fast5 directory path
  --output OUTPUT, -o OUTPUT
                        output directory path
  --cuda CUDA           choose the free GPU
  --chimeric            If provided, will output cross and within PAF
  --model {NA12878,microbe}
                        choose the pretrain model, NA12878 or microbe
```
