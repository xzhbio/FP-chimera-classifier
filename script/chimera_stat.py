import os,argparse

def reads_overlap(input):
    '''计算query之间的overlap和map长度之间的比值'''
    read_id = ''
    ratio = []
    with open(input,'r') as file:
        for reads in file:
            read = reads.split()
            if read_id == read[0]:
                if int(read[2]) > start and int(read[2]) < end:
                    overlap = end - int(read[2])
                    map_length = min(int(read[10]),length)
                    num = overlap / map_length
                    ratio.append(num)
                elif int(read[3]) > start and int(read[3]) < end:
                    overlap = int(read[3]) - start
                    map_length = min(int(read[10]),length)
                    num = overlap / map_length
                    ratio.append(num)
            else:
                read_id = read[0]
                start = int(read[2])
                end = int(read[3])
                length = int(read[10])
    return ratio

def reads_distance(input, output):
    '''计算query之间的overlap长度(碱基个数)'''
    reads_id = ''
    with open(input,'r') as f:
        lines = f.readlines()
        for line in lines:
            reads = line.split()
            if reads_id == reads[0]:
                if int(reads[2]) > start and int(reads[2]) < end:
                    distance = int(reads[2]) - end
                    with open(output,'a') as a:
                        a.write(f'{distance}\n')
                elif int(reads[3]) > start and int(reads[3]) < end:
                    distance = start - int(reads[3])
                    with open(output,'a') as a:
                        a.write(f'{distance}\n')
                else:
                    distance = min(abs(int(reads[2]) - end),abs(int(reads[3]) - start))
                    with open(output,'a') as a:
                        a.write(f'{distance}\n')
				# if (int(reads[2]) > start and int(reads[2]) < end) or (int(reads[3]) > start and int(reads[3]) < end):
				#     distance = 0
				#     with open('/home/xzh/data/chimera_anlysis/result/within_result/1_distance_reads.txt','a') as a:
				#         a.write(f'{distance}\n')
				# else:
				#     distance = min(abs(int(reads[2]) - end),abs(int(reads[3]) - start))
				#     with open('/home/xzh/data/chimera_anlysis/result/within_result/1_distance_reads.txt','a') as a:
				#         a.write(f'{distance}\n')
            else:
                reads_id = reads[0]
                start = int(reads[2])
                end = int(reads[3])
				
def ref_distance(input):
    '''计算within_chimeric mapping位点之间的distance'''
    reads_id = ''
    str_name = ''
    distance = []
    if not os.path.getsize(input):
        with open(input,'r') as f:
            lines = f.readlines()
            for line in lines:
                reads = line.split()
                if reads_id == reads[0] and str_name == reads[5]:
                    distance.append(min(abs(int(reads[7]) - end),abs(int(reads[8]) - start)))
                    # with open(output + '/within_ref_distance.txt','a') as a:
                    #     a.write(f'{distance}\n')
                else:
                    reads_id = reads[0]
                    str_name = reads[5]
                    start = int(reads[7])
                    end = int(reads[8])
    return distance

def within_strand(input,output):
    reads_id = ''
    str_name = ''
    ff = 0
    fm = 0 
    mm = 0
    mf = 0

    with open(input,'r') as f:
        lines = f.readlines()
        for line in lines:
            reads = line.split()
            if reads_id == reads[0]:
                if strand == '+' and strand == reads[4] and t == 0:
                    ff += 1
                    t = 1
                elif strand == '-'  and strand == reads[4] and a == 0:
                    mm += 1
                    a = 1
                elif g == 0:
                    fm += 1
                    g = 1
            else:
                reads_id = reads[0]
                strand = reads[4]
                str_name = reads[5]
                start = int(reads[7])
                end = int(reads[8])
                temp = reads
                t = 0
                a = 0
                g = 0

    with open(output,'w') as w:
        w.write(f'+_+:{ff}\n-_-:{mm}\n+_-:{fm}')

def hairpin_chimera(input,output):
    d = {}
    readsid = []
    with open(input,'r') as file:
        for line in file:
            reads = line.split()
            if int(reads[10]) < 500:
                continue
            id = reads[0] + "$" + reads[5]
            query_start = int(reads[2])
            query_end = int(reads[3])
            map_length = int(reads[10])
            query_length = int(reads[1])
            ref_start = int(reads[7])
            ref_end = int(reads[8])
            label = reads[4]
            if id not in d:
                d[id]={label:[[query_start,query_end,map_length,query_length,ref_start,ref_end]]}
            elif label not in d[id]:
                temp = d[id]
                temp[label] = [[query_start,query_end,map_length,query_length,ref_start,ref_end]]
            else:
                d[id][label].append([query_start,query_end,map_length,query_length,ref_start,ref_end])

    for id in d:
        if '+' in d[id] and '-' in d[id]:
            sum_a = sum(item[2] for item in d[id]['+'])
            sum_b = sum(item[2] for item in d[id]['-'])
            q_length = d[id]['+'][0][3]
            if (sum_a + sum_b) / q_length < 0.8:
                continue
            l1 = d[id]['+']
            l2 = d[id]['-']
            tag = 0
            for i in range(0,len(l1)):
                if tag == 1:
                    break
                for j in range(0,len(l2)):
                    if tag == 1:
                        break
                    # if int(l1[i][0]) < int(l2[j][0]) and int(l1[i][1]) > int(l2[j][0]) and (int(l1[i][1]) - int(l2[j][0])) > 50:
                    #     tag = 1
                    #     readsid.append(id.split('$')[0])
                    #     # with open(output + '/hairpin.txt','a') as w:
                    #     #     w.write(id.split('$')[0] + '\n')
                    #     break
                    # elif int(l1[i][0]) > int(l2[j][0]) and int(l1[i][0]) < int(l2[j][1]) and (int(l2[j][1]) - int(l1[i][0])) > 50:
                    #     tag = 1
                    #     readsid.append(id.split('$')[0])
                    #     # with open(output + '/hairpin.txt','a') as w:
                    #     #     w.write(id.split('$')[0] + '\n')
                    #     break
                    else:
                        if l1[i][0] > l2[j][0] and l1[i][0] < l2[j][1]:
                            overlap = l2[j][1] - l1[i][0]
                            overlap_ratio = overlap/min(l1[i][2],l2[j][2])
                            if overlap_ratio > 0.05:
                                continue
                        elif l1[i][1] > l2[j][0] and l1[i][1] < l2[j][1]:
                            overlap = l1[i][1] - l2[j][0]
                            overlap_ratio = overlap/min(l1[i][2],l2[j][2])
                            if overlap_ratio > 0.05:
                                continue
                        
                        # 判断ref之间的overlap_ratio
                        if l1[i][4] > l2[j][5] or l2[j][4] > l1[i][5]:
                            continue
                        if l1[i][4] > l2[j][4] and l1[i][4] < l2[j][5]:
                            overlap = l2[j][5] - l1[i][4]
                            overlap_ratio = overlap/min(l1[i][5] - l1[i][4],l2[j][5] - l2[j][4])
                            if overlap_ratio < 0.8:
                                continue
                        elif l1[i][5] > l2[j][4] and l1[i][5] < l2[j][5]:
                            overlap = l1[i][5] - l2[j][4]
                            overlap_ratio = overlap/min(l1[i][5] - l1[i][4],l2[j][5] - l2[j][4])
                            if overlap_ratio < 0.8:
                                continue
                        tag = 1
                        readsid.append(id.split('$')[0])
                        
    with open(input, 'r') as f_in, open(output, 'w') as f_out:
        for line in f_in:
            reads = line.split('\t')[0]
            if reads in readsid:
                f_out.write(line)
                    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input','-i',type=str,help='输入fastq文件路径',required=True)
    parser.add_argument('--output','-o',type=str,help='输出路径',default='.')
    args = parser.parse_args()
    hairpin_chimera(args.input, args.output)

            
            