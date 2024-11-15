# from fnmatch import fnmatch
import os
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import argparse

def calculate_base_coverage(positions):
    positions.sort(key=lambda x: x[0])  # Sort positions by start point
    covered_bases = 0
    prev_end = -1
    
    for start, end in positions:
        if start > prev_end:  # No overlap with previous interval
            covered_bases += end - start + 1
        elif (end - prev_end) > 0:  # Overlap with previous interval
            covered_bases += end - prev_end
        
        prev_end = max(prev_end, end)
    
    return covered_bases

def generate_paf(input,output,reference):
    ref_name = ''

    ref = SeqIO.to_dict(SeqIO.parse(reference,"fasta"))

    with open(input,'r') as file:
        lines = file.readlines()
        d = {}
        for line in lines:
            reads = line.split()
            if reads[0] not in d:
                d[reads[0]] = []
            d[reads[0]].append(reads)
            
    within = set()
    cross = set()

    blast_tag = 0
    for id in d:
        record = d[id]
        mapping_length = 0    
        for i in range(len(record)):
            mapping_length += int(record[i][10])
            if int(record[i][1]) == 0:
                print(record[i]) 
            if mapping_length/int(record[i][1]) < 0.8:
                continue
            read1 = record[i]
            ref_name = read1[5]
            ref_start = int(read1[7])
            ref_end = int(read1[8])
            query_start = int(read1[2])
            query_end = int(read1[3])
            length = int(read1[10])
            strand = read1[4]

            seq1 = ref[ref_name].seq[ref_start:ref_end]
            record1 = SeqRecord(seq1, id=read1[0], description="Sequence 1")
            

            for j in range(i+1, len(record)):
                reads = record[j]
                if ref_name == reads[5]:
                    if int(reads[10]) < 500:
                        continue
                    if int(reads[7]) > ref_end or int(reads[8]) < ref_start:
                        if int(reads[2]) > query_start and int(reads[2]) < query_end:
                            overlap = query_end - int(reads[2])
                            map_length = min(int(reads[10]),length)
                            num = overlap / map_length
                            if num > 0.05:
                                continue
                        elif int(reads[3]) > query_start and int(reads[3]) < query_end:
                            overlap = int(reads[3]) - query_start
                            map_length = min(int(reads[10]),length)
                            num = overlap / map_length
                            if num > 0.05:
                                continue
                        if min(abs(int(reads[7]) - ref_end),abs(int(reads[8]) - ref_start)) < 5000:
                            continue

                        seq2 = ref[reads[5]].seq[int(reads[7]):int(reads[8])]
                        record2 = SeqRecord(seq2, id="query", description="Sequence 2")
                        SeqIO.write(record1, "seq1.fasta", "fasta")
                        subprocess.run(["makeblastdb -in seq1.fasta -dbtype nucl -out seqdb > blast.log"],shell=True,check=True)
                        SeqIO.write(record2, "seq2.fasta", "fasta")
                        subprocess.run(["blastn -query seq2.fasta -db seqdb -outfmt 6 > blast_results.txt"],shell=True,check=True)
                        blast_tag = 1
                        blast_query = []
                        blast_ref = []
                        if os.path.getsize('blast_results.txt'):
                            with open('blast_results.txt','r') as f:
                                for row in f:
                                    blast_record = row.split()
                                    name = blast_record[1]
                                    length = int(blast_record[3])
                                    q_start = int(blast_record[6])
                                    q_end = int(blast_record[7])
                                    r_start = int(blast_record[8])
                                    r_end = int(blast_record[9])
                                    if r_start > r_end:
                                        t = r_start
                                        r_start = r_end
                                        r_end = t

                                    blast_query.append([q_start,q_end])
                                    blast_ref.append([r_start,r_end])
                            len_q = calculate_base_coverage(blast_query)
                            len_r = calculate_base_coverage(blast_ref)
                            map_length = min(len_q,len_r)
                            ratio = map_length / min(length,int(reads[10]))
                            if ratio > 0.05:
                                continue
                        if strand == reads[4]:
                            within.add('\t'.join(reads))
                            within.add('\t'.join(read1))
                        else:
                            within.add('\t'.join(reads))
                            within.add('\t'.join(read1))
                        # with open(output + '/within.paf','a') as w1:
                        #     w1.write('\t'.join(reads) + '\n')
                        #     if within==0:
                        #         w1.write('\t'.join(read1) + '\n')
                        #         within = 1


                elif ref_name != reads[5]:
                    if int(reads[10]) < 500:
                        continue
                    if int(reads[2]) > query_start and int(reads[2]) < query_end:
                        overlap = query_end - int(reads[2])
                        map_length = min(int(reads[10]),length)
                        num = overlap / map_length
                        if num > 0.05:
                            continue
                    elif int(reads[3]) > query_start and int(reads[3]) < query_end:
                        overlap = int(reads[3]) - query_start
                        map_length = min(int(reads[10]),length)
                        num = overlap / map_length
                        if num > 0.05:
                            continue
                    seq2 = ref[reads[5]].seq[int(reads[7]):int(reads[8])]
                    record2 = SeqRecord(seq2, id="query", description="Sequence 2")
                    SeqIO.write(record1, "seq1.fasta", "fasta")
                    subprocess.run(["makeblastdb -in seq1.fasta -dbtype nucl -out seqdb > blast.log"],shell=True,check=True)
                    SeqIO.write(record2, "seq2.fasta", "fasta")
                    subprocess.run(["blastn -query seq2.fasta -db seqdb -outfmt 6 > blast_results.txt"],shell=True,check=True)
                    blast_tag = 1
                    blast_query = []
                    blast_ref = []
                    if os.path.getsize('blast_results.txt'):
                        with open('blast_results.txt','r') as f:
                            for row in f:
                                blast_record = row.split()
                                name = blast_record[1]
                                length = int(blast_record[3])
                                q_start = int(blast_record[6])
                                q_end = int(blast_record[7])
                                r_start = int(blast_record[8])
                                r_end = int(blast_record[9])
                                if r_start > r_end:
                                    t = r_start
                                    r_start = r_end
                                    r_end = t

                                blast_query.append([q_start,q_end])
                                blast_ref.append([r_start,r_end])
                        len_q = calculate_base_coverage(blast_query)
                        len_r = calculate_base_coverage(blast_ref)
                        map_length = min(len_q,len_r)
                        ratio = map_length / min(length,int(reads[10]))
                        if ratio > 0.05:
                            continue

                    if strand == reads[4]:
                        cross.add('\t'.join(reads))
                        cross.add('\t'.join(read1))
                    else:
                        cross.add('\t'.join(reads))
                        cross.add('\t'.join(read1))
                    # with open(output + '/cross.paf','a') as w2:
                    #     w2.write('\t'.join(reads) + '\n')
                    #     if cross==0:
                    #         w2.write('\t'.join(read1) + '\n')
                    #         cross = 1
    if blast_tag == 1:
        subprocess.run(['rm seqdb.* seq1.fasta seq2.fasta blast_results.txt'],shell=True,check=True)

    with open(output + '/within.paf','w') as w:
        within_str = '\n'.join(within) + '\n'
        w.write(within_str)

    with open(output + '/cross.paf','w') as w:
        cross_str = '\n'.join(cross) + '\n'
        w.write(cross_str)
       
    # print(os.path.join(output,'/within.paf'))
    os.system(f"sort {os.path.join(output,'within.paf')} -o {os.path.join(output,'within.paf')}")
    os.system(f"sort {os.path.join(output,'cross.paf')} -o {os.path.join(output,'cross.paf')}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input','-i',type=str,help='input path of fastq file',required=True)
    parser.add_argument('--ref',type=str,help='reference path',required=True)
    parser.add_argument('--output','-o',type=str,help='output directory path',default='.')
    args = parser.parse_args()

    generate_paf(args.input, args.output, args.ref)