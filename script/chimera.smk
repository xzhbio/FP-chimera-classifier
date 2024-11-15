import os
out_dir = config['output']
dir_path = config['dir_path']
if config['model'] == 'NA12878':
    model = os.path.join(dir_path,'pretrain_model/ResNet_NA12878.pth')
else:
    model = os.path.join(dir_path,'pretrain_model/ResNet_micro_all.pth')

rule all:
    input:
        expand("{out_dir}/align.paf",out_dir=out_dir),
        expand("{out_dir}/FP-SVs_reads_id.txt",out_dir=out_dir),
        (expand("{out_dir}/cross.paf",out_dir=out_dir) if config['chimeric'] else [] +
        expand("{out_dir}/within.paf",out_dir=out_dir) if config['chimeric'] else [])

rule minimap2:
    input:
        fastq = config['input'],
        ref = config['ref']
    output:
        "{out_dir}/align.paf"
    shell:
        "minimap2 -x map-ont -t 24 -c {input.ref} {input.fastq} | sort -k 1,1 -k 6,6 | awk '$12>=60' > {output}"

rule chimeric:
    input:
        "{out_dir}/align.paf"
    params:
        config['ref']
    output:
        "{out_dir}/cross.paf",
        "{out_dir}/within.paf"
    shell:
        "python {dir_path}/script/chimera_paf.py -i {input} --ref {params} --output {out_dir}"

rule FPSV_detect:
    input:
        "{out_dir}/align.paf"
    output:
        "{out_dir}/hairpin.paf"
    shell:
        "python {dir_path}/script/chimera_stat.py --input {input} --output {output}"

rule signal_data_prepare:
    input:
        fast5 = config['fast5_path'],
        hairpin = "{out_dir}/hairpin.paf"
    output:
        "{out_dir}/hairpin.npy"
    shell:
        "python {dir_path}/script/data_prepare.py --fast5 {input.fast5} --input_path {input.hairpin} --save {output} --length 20000"

rule model_check:
    input:
        external = "{out_dir}/hairpin.npy"
    output:
        "{out_dir}/FP-SVs_reads_id.txt"
    params:
        cuda = config['cuda']
    log:
        "{out_dir}/model_check.log"
    # conda:
    #     "pytorch_cu11.6.yaml"
    shell:
        "python {dir_path}/script/signal_classify_ResNet.py --load {model} --external {input.external} --output {output} --cuda {params.cuda}"