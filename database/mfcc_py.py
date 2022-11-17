import librosa
import os
import numpy as np
np.set_printoptions(threshold = np.inf)
import csv

from cProfile import label
from ctypes import sizeof

##MFCC feature extract
##
database_path = r"G:\py_mfcc_n_dtw\database\\"
imfor_dict = [
    #gender: 0 stands for male, 1 stands for female
    {'gender':0, 'age':30, 'accent':'german'},
    {'gender':0, 'age':25, 'accent':'german'},
    {'gender':0, 'age':31, 'accent':'german'},
    {'gender':0, 'age':23, 'accent':'german'},
    {'gender':0, 'age':25, 'accent':'german'},
    {'gender':0, 'age':25, 'accent':'german'},
    {'gender':0, 'age':27, 'accent':'german/spanish'},
    {'gender':0, 'age':41, 'accent':'german'},
    {'gender':0, 'age':35, 'accent':'south korea'},
    {'gender':0, 'age':36, 'accent':'german'},
    {'gender':0, 'age':33, 'accent':'german'},
    {'gender':1, 'age':26, 'accent':'german'},
    {'gender':0, 'age':27, 'accent':'german'},
    {'gender':0, 'age':31, 'accent':'spanish'},
    {'gender':0, 'age':28, 'accent':'madras'},
    {'gender':0, 'age':30, 'accent':'german'},
    {'gender':0, 'age':26, 'accent':'german'},
    {'gender':0, 'age':25, 'accent':'levant'},
    {'gender':0, 'age':23, 'accent':'english'},
    {'gender':0, 'age':25, 'accent':'german'},
    {'gender':0, 'age':26, 'accent':'german'},
    {'gender':0, 'age':33, 'accent':'german'},
    {'gender':0, 'age':28, 'accent':'german'},
    {'gender':0, 'age':26, 'accent':'chinese'},
    {'gender':0, 'age':22, 'accent':'brasilian'},
    {'gender':1, 'age':22, 'accent':'chinese'},
    {'gender':0, 'age':31, 'accent':'italian'},
    {'gender':1, 'age':28, 'accent':'german'},
    {'gender':0, 'age':23, 'accent':'german'},
    {'gender':0, 'age':28, 'accent':'german'},
    {'gender':0, 'age':26, 'accent':'german'},
    {'gender':0, 'age':23, 'accent':'egyptian_american'},
    {'gender':0, 'age':26, 'accent':'german'},
    {'gender':0, 'age':25, 'accent':'german'},
    {'gender':0, 'age':24, 'accent':'chinese'},
    {'gender':0, 'age':22, 'accent':'german'},
    {'gender':0, 'age':27, 'accent':'italian'},
    {'gender':0, 'age':32, 'accent':'spanish'},
    {'gender':0, 'age':29, 'accent':'german'},
    {'gender':0, 'age':26, 'accent':'german'},
    {'gender':0, 'age':30, 'accent':'south african'},
    {'gender':0, 'age':29, 'accent':'arabic'},
    {'gender':1, 'age':31, 'accent':'german'},
    {'gender':0, 'age':61, 'accent':'german'},
    {'gender':0, 'age':'nan', 'accent':'german'},
    {'gender':0, 'age':30, 'accent':'german'},
    {'gender':1, 'age':23, 'accent':'danish'},
    {'gender':0, 'age':26, 'accent':'german'},
    {'gender':0, 'age':26, 'accent':'german'},
    {'gender':0, 'age':24, 'accent':'german'},
    {'gender':0, 'age':26, 'accent':'german'},
    {'gender':1, 'age':34, 'accent':'french'},
    {'gender':0, 'age':24, 'accent':'german'},
    {'gender':0, 'age':27, 'accent':'german'},
    {'gender':0, 'age':23, 'accent':'german'},
    {'gender':1, 'age':24, 'accent':'german'},
    {'gender':1, 'age':27, 'accent':'german'},
    {'gender':1, 'age':29, 'accent':'german'},
    {'gender':1, 'age':31, 'accent':'german'},
    {'gender':1, 'age':27, 'accent':'tamil'}
]
#print(len(imfor_dict))

for ii in range(0,10):
    #find_dir  = r"G:\py_mfcc_n_dtw\database\mfcc_csv\find_directory" + r"\\" + str(ii) + ".csv"
    #query_dir = r"G:\py_mfcc_n_dtw\database\mfcc_csv\query" + r"\\" + str(ii) + ".csv"
    #find_dir  = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut\find_directory" + r"\\" + str(ii) + ".csv"
    #query_dir = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut\query" + r"\\" + str(ii) + ".csv"
    find_dir  = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut_cut\find_directory" + r"\\" + str(ii) + ".csv"
    query_dir = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut_cut\query" + r"\\" + str(ii) + ".csv"
    col_name_data = ['gender', 'age', 'accent', 'mfcc_feature']
    f1 = open(find_dir, mode = 'w', encoding ='utf-8', newline = "")
    csv_writer = csv.writer(f1)
    csv_writer.writerow(col_name_data)
    f2 = open(query_dir, mode = 'w', encoding ='utf-8', newline = "")
    csv_writer = csv.writer(f2)
    csv_writer.writerow(col_name_data)

#for j in range(1,46):
#for j in range(1,21): #cut half
for j in range(1,11): #cut three of quarter
    if(j < 10):
        browse_name = '0'+str(j)
    else:
        browse_name = str(j)
    data_path = database_path + browse_name + r"\\"
    file_names = os.listdir(data_path)
    for i in range(0, len(file_names)):
        temp_path = data_path + file_names[i]
        #Load an audio file as a floating point time series.
        ##sr: sample rates, to preserve the native sampling rate of the file, use sr=None.
        y,sr = librosa.load(temp_path, sr=None)
        ##label the number of speech audio
        label = file_names[i][0]
        ##extract MFCC feature, 13 dimension
        mfcc_feature = librosa.feature.mfcc(y=y, sr=sr, n_mfcc=13)
        #transpose for DTW
        mfcc_feature_t = np.transpose(mfcc_feature)
        ##writing csv
        #csv_path = r"G:\py_mfcc_n_dtw\database\mfcc_csv\find_directory" + r"\\" + label + ".csv"
        #csv_path = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut\find_directory" + r"\\" + label + ".csv"
        csv_path = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut_cut\find_directory" + r"\\" + label + ".csv"
        w_data = [imfor_dict[j-1]['gender'], imfor_dict[j-1]['age'], imfor_dict[j-1]['accent'], mfcc_feature_t]

        #print("labeled data is")
        #print(w_data)

        f = open(csv_path, mode = 'a', encoding ='utf-8', newline = "")
        csv_writer = csv.writer(f)

        csv_writer.writerow(w_data)

        f.close()

#for k in range(46,61):
#for k in range(21,31):
for k in range(11,16):
    if(k < 10):
        browse_name = '0'+str(k)
    else:
        browse_name = str(k)
    data_path = database_path + browse_name + r"\\"
    file_names = os.listdir(data_path)
    for m in range(0, len(file_names)):
        temp_path = data_path + file_names[m]
        #Load an audio file as a floating point time series.
        ##sr: sample rates, to preserve the native sampling rate of the file, use sr=None.
        y,sr = librosa.load(temp_path, sr=None)
        ##label the number of speech audio
        label = file_names[m][0]
        ##extract MFCC feature, 13 dimension
        mfcc_feature = librosa.feature.mfcc(y=y, sr=sr, n_mfcc=13)
        #transpose for DTW
        mfcc_feature_t = np.transpose(mfcc_feature)
        ##writing csv
        #csv_path = r"G:\py_mfcc_n_dtw\database\mfcc_csv\query" + r"\\" + label + ".csv"
        #csv_path = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut\query" + r"\\" + label + ".csv"
        csv_path = r"G:\py_mfcc_n_dtw\database\mfcc_csv_cut_cut\query" + r"\\" + label + ".csv"
        w_data = [imfor_dict[k-1]['gender'], imfor_dict[k-1]['age'], imfor_dict[k-1]['accent'], mfcc_feature_t]

        #print("labeled data is")
        #print(w_data)

        f = open(csv_path, mode = 'a', encoding ='utf-8', newline = "")
        csv_writer = csv.writer(f)

        csv_writer.writerow(w_data)

        f.close()