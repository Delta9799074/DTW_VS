import time
import numpy as np
np.set_printoptions(threshold = np.inf)
import random
import pandas as pd
from fastdtw import fastdtw
from dtaidistance import dtw_visualisation as dtwvis
from scipy.spatial.distance import euclidean

def compare_feature(library_path, number, cmp_feature):
    csv_path = library_path + r"\\" + str(number) + '.csv'
    pdreader0 = pd.read_csv(csv_path)
    pd_mfcc_data0 = pdreader0['mfcc_feature']

    ''' y,sr = librosa.load("G:\\hls\\py_mfcc_n_dtw\\database\\01\\0_01_3.wav", sr=None)
    mfcc_feature = librosa.feature.mfcc(y=y, sr=sr, n_mfcc=13)
    print(mfcc_feature.shape) '''

    wrong_cnt0 = 0
    feature0_read_correct = False
    distance = []

    for i0 in range(0, len(pd_mfcc_data0)):
        temp_str = pd_mfcc_data0[i0].strip('[')
        temp_str = temp_str.strip(']')
        temp_str = temp_str.replace('\n', '')
        temp_str = temp_str.replace('] [', '\n')
        nlines = temp_str.count('\n')
        nlines = nlines + 1
        temp_str_list = temp_str.split('\n')
        np_mfcc_feature = np.empty([nlines, 13], dtype = float)
        this_line_wrong_tag = 0

        for j in range(0,nlines):
            line_list = temp_str_list[j].split(' ')
            value_moved = ''
            line_list_filted = list(filter(lambda val: val != value_moved, line_list))
            if(len(line_list_filted) != 13):
                this_line_wrong_tag += 1
            else:
                for ii in range(0,13):
                    np_mfcc_feature[j][ii] = eval(line_list_filted[ii])
        if(this_line_wrong_tag != 0):
            wrong_cnt0 += 1
        else:
            temp_distance, temp_path = fastdtw(cmp_feature, np_mfcc_feature, dist = euclidean)
            distance.append(temp_distance)

    if(wrong_cnt0 == 0):
        feature0_read_correct = True
    
    if(feature0_read_correct):
        result = min(distance)
    
    return result

def gen_query_feature(qlibrary_path, number):
    query_path = qlibrary_path + r"\\" + str(number) + '.csv'
    pdreader_query0 = pd.read_csv(query_path)
    pd_mfcc_query_data0 = pdreader_query0['mfcc_feature']
    feature0_read_correct = False
    query_index = random.randint(0,(len(pd_mfcc_query_data0)-1))

    temp_query_str = pd_mfcc_query_data0[query_index].strip('[')
    temp_query_str = temp_query_str.strip(']')
    temp_query_str = temp_query_str.replace('\n', '')
    temp_query_str = temp_query_str.replace('] [', '\n')
    nlines = temp_query_str.count('\n')
    nlines = nlines + 1
    temp_query_str_list = temp_query_str.split('\n')
    np_mfcc_query_feature = np.empty([nlines, 13], dtype = float)
    this_line_wrong_tag = 0
    for j in range(0,nlines):
        line_list = temp_query_str_list[j].split(' ')
        value_moved = ''
        line_list_filted = list(filter(lambda val: val != value_moved, line_list))
        if(len(line_list_filted) != 13):
            this_line_wrong_tag += 1
        else:
            for ii in range(0,13):
                np_mfcc_query_feature[j][ii] = eval(line_list_filted[ii])
    return np_mfcc_query_feature

def query_func(lib_path, np_mfcc_query_feature):
    query_distance =[]
    query_time = []
    for q in range(0,10):
        query_time_start = time.time()
        distance_query = compare_feature(lib_path, q, np_mfcc_query_feature)
        query_time_end = time.time()
        query_time_temp = query_time_end - query_time_start
        #print('query', str(q), 'execute time:', query_time_temp)
        query_distance.append(distance_query)
        query_time.append(query_time_temp)
    min_dist = min(query_distance)
    min_idx = query_distance.index(min_dist)

    np_exe_time = np.array(query_time)
    return min_idx,np_exe_time

qlib_path = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut\query"
lib_path  = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut\find_directory"

#query_path = r"G:\py_mfcc_n_dtw\database\mfcc_csv\query\\0.csv"
#pdreader_query0 = pd.read_csv(query_path) #gen length

correct_num = 0
wrong_num = 0

tt_time = []

for i in range(0,1):
#for i in range(0,pdreader_query0.shape[0]):
    query_feat = gen_query_feature(qlib_path,0)
    query_idx, query_time = query_func(lib_path, query_feat)
    if(query_idx == 0):
        correct_num += 1
    else:
        wrong_num += 1
    tt_time.append(query_time)

np_tt_time = np.array(tt_time)
correct_percent = correct_num / (correct_num + wrong_num)

print("RECOGNITION CORRECT RATE:", correct_percent)
print("TOTAL EXECUTE TIME:", np_tt_time.sum())