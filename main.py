import argparse, subprocess
# import snakemake
from script import chimera_stat as stat,chimera_paf
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='FP-SVs classifier')
    parser.add_argument('--input','-i',type=str,help='input path of fastq file',required=True)
    parser.add_argument('--ref',type=str,help='reference path',required=True)
    parser.add_argument('--fast5',type=str,help='fast5 directory path',required=True)
    parser.add_argument('--output','-o',type=str,help='output directory path',default='.')
    parser.add_argument('--cuda',type=str,help='choose the free GPU',default='0')
    parser.add_argument('--chimeric',action='store_true',help='If provided, will output cross and within PAF')
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    parser.add_argument('--model',choices=['NA12878','microbe'],help='choose the pretrain model, NA12878 or microbe')
    args = parser.parse_args()
    # _,file_extension = os.path.splitext(args.input)
    # if file_extension != '.paf':
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    if args.chimeric:
        chimeric=1
    else:
        chimeric=0

    command = ['snakemake', '--snakefile', f'{os.path.join(current_dir,"script/chimera.smk")}', '--config', f'input={args.input}', f'dir_path={current_dir}',
               f'ref={args.ref}', f'output={args.output}', f'fast5_path={args.fast5}', f'cuda={args.cuda}', f'model={args.model}', f'chimeric={chimeric}', '--cores', 'all']
    
    subprocess.run(command)