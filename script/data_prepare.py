import os
import h5py
import numpy as np
import random
import argparse

def outside_filt(input_list):
    x = np.where(np.array(input_list) > 3*np.median(input_list))
    x = x[0]
    if len(x) > 0:
        for i in range(len(x)):
            # print(x[i])
            if x[i] == 0:
                input_list[x[i]] = input_list[x[i] + 1]
            elif x[i] == len(input_list) - 1:
                input_list[x[i]] = input_list[x[i] - 1]
            else:
                input_list[x[i]] = (input_list[x[i] - 1] + input_list[x[i] + 1]) / 2
    return input_list

def generate_chimeric_breakpoint(chimeric_path_list):
    # 返回嵌合序列id及其嵌合断点以生成电信号数据集
    chimeric_breakpoint_dict = {}
    for chimeric_path in chimeric_path_list:
        tag = 0
        with open(chimeric_path, 'r') as f:
            for row in f:
                read = row.split()
                if tag == 0:
                    query_id = read[0]
                    query_start = int(read[2])
                    tag = 1
                if query_id not in chimeric_breakpoint_dict:
                    chimeric_breakpoint_dict[query_id] = []
                if read[0] != query_id:
                    query_id = read[0]
                    query_start = int(read[2])
                elif int(read[2]) == query_start:
                    continue
                elif int(read[2]) > query_start:
                    chimeric_breakpoint_dict[query_id].append(int(read[2]))
                else:
                    chimeric_breakpoint_dict[query_id].append(int(query_start))
    return chimeric_breakpoint_dict

def traverse_directory(directory):
    file_paths = []
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if file_path.endswith('.fast5'):
                file_paths.append(file_path)
    return file_paths

def standardization(input):
    mean_signal = np.mean(input)
    std_signal = np.std(input)

    if std_signal == 0:
        signal = input - mean_signal
    else:
        signal = (input - mean_signal) / std_signal
    return signal

def generate_signal_dataset(fast5_path, chimeric_path_list, output, length):
    chimeric_breakpoint_dict = generate_chimeric_breakpoint(chimeric_path_list)
    fast5_file_list = traverse_directory(fast5_path)
    signal_data = []
    filt_num = 0
    pad = int(length / 2)
    for fast5_flie in fast5_file_list:
        try:
            fast5 = h5py.File(fast5_flie)
            for reads in fast5.keys():
                reads_id = reads.split('_')[1]
                if reads_id in chimeric_breakpoint_dict:
                    for breakpoint in chimeric_breakpoint_dict[reads_id]:
                        # if breakpoint*10 + 10000 > len(fast5[f'{reads}']['Raw']['Signal'][:]) or breakpoint*10 - 10000 < 0:
                        #     print('fail')
                        #     continue
                        if breakpoint*10 - pad > len(fast5[f'{reads}']['Raw']['Signal'][:]):
                            filt_num += 1
                            continue
                        if breakpoint*10 - pad > 0:
                            signal = np.array(fast5[f'{reads}']['Raw']['Signal'][breakpoint*10 - pad: breakpoint*10 + pad])
                        else:
                            signal = np.array(fast5[f'{reads}']['Raw']['Signal'][0: breakpoint*10 + pad])
                        # 离群值过滤
                        signal = outside_filt(signal)
                        # z-score 标准化
                        signal = standardization(signal)
                        if len(signal) < 2*pad:
                            signal = np.pad(signal, (0, 2*pad - len(signal)), mode='constant', constant_values=0.0)
                        label = 1
                        signal_data.append([signal, label, reads_id])
        except Exception as e:
            print(fast5_flie,e)
    signal_data = np.array(signal_data, dtype=object)
    np.save(output, signal_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fast5',type=str,required=True)
    parser.add_argument('--input_path','-i',type=str,required=True)  # path of hairpin.paf
    parser.add_argument('--save',type=str)
    parser.add_argument('--length',type=int)
    args = parser.parse_args()
    generate_signal_dataset(args.fast5, [args.input_path], args.save, args.length)
 