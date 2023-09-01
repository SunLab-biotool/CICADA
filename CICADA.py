###################################################################### -*-coding : utf-8-*-
import re
import os
import joblib
import pandas as pd
import optparse
from optparse import OptionParser
import time
import subprocess
import math
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import sys
import faulthandler;faulthandler.enable()
import uuid
import shutil

########################################################################
def rfPredict(data):
    rf = joblib.load(curPath + '/CICADA_Parameters/circRNA_RF_Predict.model')
    pre = rf.predict_proba(data)
    return pre

########################################################################
def translate(dna):
    gencode = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W'}
    amino_acid_sequence = ""
    for start in range(0, len(dna) - 2, 3):
        stop = start + 3
        codon = dna[start:stop]
        aa = gencode.get(codon.upper(), 'X')
        amino_acid_sequence = amino_acid_sequence + aa
    return amino_acid_sequence

########################################################################
def write_table_colname(table_name,label):
    if label == 1:
        custom_cols = ["circRNA ID","Strand","circRNA sequence","circRNA length","circHORF sequence","circHORF score","circHORF start","circHORF end","circHORF length","Product sequence"]
    else:
        custom_cols = ["circRNA ID","Coding probability score","HPCR sequence","HPCR length","HPCR score","HPCR coverage","Fickett score","Conservation","m6A number"]
    table_name.loc[-1] = custom_cols
    table_name.index = table_name.index + 1
    table_name.sort_index(inplace=True)
    return table_name
########################################################################
def get_maohao_info(table_name):
    table_name.iloc[:,0] = table_name.iloc[:,0].str.replace('.*:', '')
    table_name.iloc[:,1] = table_name.iloc[:,1].str.replace('.*:', '')
    table_name.iloc[:,2] = table_name.iloc[:,2].str.replace('.*:', '')
    table_name.iloc[:,3] = table_name.iloc[:,3].str.replace('.*:', '')
    table_name.iloc[:,4] = table_name.iloc[:,4].str.replace('.*:', '')
    table_name.iloc[:,5] = table_name.iloc[:,5].str.replace('.*:', '')
    table_name.iloc[:,6] = table_name.iloc[:,6].str.replace('.*:', '')
    table_name.iloc[:,7] = table_name.iloc[:,7].str.replace('.*:', '')
    table_name.iloc[:,8] = table_name.iloc[:,8].str.replace('.*:', '')
    table_name.iloc[:,9] = table_name.iloc[:,9].str.replace('.*:', '')
#    table_name.iloc[:,10] = table_name.iloc[:,10].astype(str)
#    table_name.iloc[:,11] = table_name.iloc[:,11].astype(str)
#    table_name.iloc[:,12] = table_name.iloc[:,12].str.replace('.*:', '')
    table_name = write_table_colname(table_name,1)
    return table_name

########################################################################
def Pre_mian(Out_Feature1, Out_Feature2, Out_Feature3, Out_Information1, Out_Information2, Out_Information3, top_num):
    try:
        data3 = pd.read_csv(Out_Feature3,sep=' ',header=None)
        data3 = data3.fillna('null')
        data3[8] = 'null'
        data3[9] = 'null'
        data3[10] = 'null'
        data3.to_csv(Out_Feature3,sep='\t', index=False, header=False)
        information3 = pd.read_csv(Out_Information3,sep='\t',header=None)
        information3[13] = 'null'
        information3[14] = 'null'
        information3[16] = 'amino_acid_sequence:null'
    except pd.errors.EmptyDataError:
        data3 = pd.DataFrame()
        information3 = pd.DataFrame()
    try:
        data2 = pd.read_csv(Out_Feature2,sep=' ',header=None)
        data2 = data2.fillna('null')
        data2[8] = 'null'
        data2[9] = 'null'
        data2[10] = [str(translate(str(i))) for i in data2[7]]
        data2.to_csv(Out_Feature2,sep='\t', index=False, header=False)
        information2 = pd.read_csv(Out_Information2,sep='\t',header=None)
        information2[13] = 'null'
        information2[14] = 'null'
        information2[16] = 'amino_acid_sequence:null'
    except pd.errors.EmptyDataError:
        data2 = pd.DataFrame()
        information2 = pd.DataFrame()
    try:
        data1 = pd.read_csv(Out_Feature1,sep=' ',header=None)
        dataused = data1.iloc[:,[1,2,3,4,5,6]].values
        pre = rfPredict(dataused)
        data1[8] = pd.Series(pre[:,0])
        data1[9] = pd.Series(pre[:,1])
        data1[10] = [str(translate(str(i))) for i in data1[7]]
        data1.to_csv(Out_Feature1,sep='\t', index=False, header=False)
        id_dict = dict(map(lambda x:[x[0].split('|')[0],(x[8],x[9])],data1.values.tolist()))
        information1 = pd.read_csv(Out_Information1,sep='\t',header=None)
        information1[13] = information1[0].apply(lambda x:id_dict.get(x.split(':')[1])[0])
        information1[14] = information1[0].apply(lambda x:id_dict.get(x.split(':')[1])[1])
        information1[16] = ['amino_acid_sequence:'+str(translate(str(i.split(':')[1]))) for i in information1[3]]
        re1_1 = information1[information1.iloc[:, 2].str.split(':').str[-1] == "-"].copy()
        re1_2 = information1[information1.iloc[:, 2].str.split(':').str[-1] == "+"].copy()
        origin_Seq_len = re1_1.iloc[:, 10].str.split(':').str[-1].astype(int)
        sequence_fold = re1_1.iloc[:, 4].str.split(':').str[-1].astype(int)
        re1_1.iloc[:, 6] = origin_Seq_len * sequence_fold - re1_1.iloc[:, 6].str.split(':').str[-1].astype(int) + 1
        re1_1.iloc[:, 6] = 'Start:' + re1_1.iloc[:, 6].astype(str)
        re1_1.iloc[:, 7] = origin_Seq_len * sequence_fold - re1_1.iloc[:, 7].str.split(':').str[-1].astype(int) + 1
        re1_1.iloc[:, 7] = 'End:' + re1_1.iloc[:, 7].astype(str)
        re2 = pd.concat([re1_2, re1_1], ignore_index=True)
        for i in range(re2.shape[0]):
            if int(re2.loc[i, 11].split(':')[-1]) < 3:
                s1 = int(re2.loc[i, 6].split(':')[-1])
                e1 = int(re2.loc[i, 7].split(':')[-1])
            else:
                e1 = int(re2.loc[i, 6].split(':')[-1])
                s1 = int(re2.loc[i, 7].split(':')[-1])
            length = int(re2.loc[i, 10].split(':')[-1])
            re2.loc[i, 6] = 'Start:' + str(int(s1 if s1 <= length else s1 - math.ceil((s1 - length) / length) * length))
            re2.loc[i, 7] = 'End:' + str(int(e1 if e1 <= length else e1 - math.ceil((e1 - length) / length) * length))
        information1 = re2
    except pd.errors.EmptyDataError:
        data1 = pd.DataFrame()
        information1 = pd.DataFrame()
    data = pd.concat([data1, data2, data3], ignore_index=True)
    data = data.iloc[:,[0,9,7,3,2,4,6,1,5]]
    data = write_table_colname(data,2)
    data.to_csv(Out_Feature,sep='\t', index=False, header=False)

    information = pd.concat([information1, information2, information3], ignore_index=True)
    information = information.iloc[:,[0,2,12,10,3,5,6,7,9,15]]
    information = get_maohao_info(information)
    if top_num == 'all':
        information_top_all = information.groupby([0], sort=False).apply(lambda x: x.sort_values([5], ascending=False)).reset_index(drop=True)
        information_top_all.to_csv(Out_Information,sep='\t',index=False,header=False)
    else:
        top_num = int(top_num)
        information_top_num = information.groupby([0], sort=False).apply(lambda x: x.sort_values([5], ascending=False)).reset_index(drop=True)
        information_top = information_top_num.groupby([0]).head(top_num)
        information_top.to_csv(Out_Information,sep='\t',index=False,header=False)
    end = time.time()

#######################################################################
def cur_file_dir():
    path = sys.path[0]
    if os.path.isdir(path):
       return path
    elif os.path.isfile(path):
       return os.path.dirname(path)
#######################################################################
def TwoLineFasta (Seq_Array):
    Tmp_sequence_Arr = []
    Tmp_trans_str = ''
    for i in range(len(Seq_Array)):
        if '>' in Seq_Array[i]:
            if i == 0:
                Tmp_sequence_Arr.append(Seq_Array[i])
            else:
                Tmp_sequence_Arr.append(Tmp_trans_str)
                Tmp_sequence_Arr.append(Seq_Array[i])
                Tmp_trans_str = ''
        else:
            if i == len(Seq_Array) - 1:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
                Tmp_sequence_Arr.append(Tmp_trans_str)
            else:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
    return Tmp_sequence_Arr
#######################################################################
def Tran_checkSeq (input_arr):
    label_Arr = []
    FastA_seq_Arr = []
    for n in range(len(input_arr)):
        if n == 0 or n % 2 == 0:
            label = input_arr[n]
            label_Arr.append(label)
        else :
            seq = input_arr[n]
            FastA_seq_Arr.append(seq)
    num = 0
    for i in range(len(label_Arr)):
        Label = label_Arr[num]
        Seq = FastA_seq_Arr[num]
        tran_fir_seq = Seq.lower()
        tran_sec_seq_one = tran_fir_seq.replace('u','t')
        tran_sec_seq = tran_sec_seq_one.replace('\r','')
        if 'n' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (n),please checkout your sequence again' + '\n'
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'w' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (w),please checkout your sequence again' + '\n'
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'd' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (d),please checkout your sequence again' + '\n'
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'r' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (r),please checkout your sequence again' + '\n'
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 's' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (s),please checkout your sequence again' + '\n'
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'y' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (y),please checkout your sequence again' + '\n'
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'm' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (m),please checkout your sequence again' + '\n'
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        num = int(num) + int(1)
    return (label_Arr,FastA_seq_Arr)
#######################################################################
def InitCodonSeq(num,length,step,Arr):
    TempStrPar = ''
    for i in range(num,length,step):
        index = i
        code1 = Arr[index]
        index += 1
        code2 = Arr[index]
        index += 1
        code3 = Arr[index]
        Temp = code1+code2+code3
        TempStrPar = TempStrPar+Temp+' '
    return TempStrPar
#######################################################################
def SequenceProcessing(sequence):
    length = len(sequence)
    tran_lower_seq = sequence.lower()
    tran_cis_seq = tran_lower_seq.replace('u','t')
    cis_sequence_Arr = list(tran_cis_seq)
    return cis_sequence_Arr
#######################################################################
def Reading(read_file):
    label_Arr_tmp = []
    FastA_seq_Arr_tmp = []
    for i in range(len(read_file)):
        if i == 0 or i % 2 == 0:
            label = read_file[i]
            label_Arr_tmp.append(label)
        else:
            seq = read_file[i]
            FastA_seq_Arr_tmp.append(seq)
    return(label_Arr_tmp,FastA_seq_Arr_tmp)
#######################################################################
def CircSimulation(line_sequence,hash_score):
    CodonScore = []
    TempArray = line_sequence.split(' ')
    TempArray.pop()
    seqLength = len(TempArray)
    WindowStep = 20
    WinLen = seqLength - WindowStep
    WinLen = WinLen +1
    for j in range(WinLen):
        number = 0
        window_seq = []
        for n in range(j,WindowStep+j):
            window_seq.append(TempArray[n])
        for t in range(0,len(window_seq)-1):
            temp1 = window_seq[t] + window_seq[t+1]
            num_temp = re.compile('[atcg]{6}')
            if num_temp.match(temp1):
                number = float(number) + float(hash_score[temp1])
        number = number / WindowStep
        CodonScore.append(number)
    return CodonScore
######################################################################
def MaxSubseqPlus(string_score,Frame_string):
    TempArray = Frame_string.split(' ')
    TempArray.pop()
    Total_Start = 0
    Total_End = 0
    Max = 0
    frame_max_string = ''
    for i in range(len(string_score)):
        sumNum = 0
        for j in range(i,len(string_score)):
            sumNum = sumNum + float(string_score[j])
            if sumNum > Max:
                Total_Start = i
                Total_End = j
                Max = sumNum
    for n in range(Total_Start, Total_End + 20):
        frame_max_string = frame_max_string + TempArray[n] + ' '
    return (Max, frame_max_string, Total_Start, Total_End)
######################################################################
def FindStartCodne(TT_seq,FT_seq,ST_seq,ET_seq,Start,End,j,Seq_len):
    TT_seq_Array = TT_seq.split(' ')
    TT_seq_Array.pop()
    ST_seq_Array = ST_seq.split(' ')
    ST_seq_Array.pop()
    FT_seq_Array = FT_seq.split(' ')
    FT_seq_Array.pop()
    ET_seq_Array = ET_seq.split(' ')
    ET_seq_Array.pop()
    Start_coden_Array = []
    End_coden_Array_left = []
    End_coden_Array_right = []
    cds_end = End * 3 + 20 * 3 + j
    if Seq_len % 3 == 0 and cds_end <= Seq_len:
        seq_fold = 3
        begin_one = int(Start + Seq_len/3)
        end_one = int(Seq_len/3 * 2 + Start + 1)
        for i in range(begin_one,end_one):
            end_temp_right1 = TT_seq_Array[i]
            if str(end_temp_right1) == 'tag' or str(end_temp_right1) == 'taa' or str(end_temp_right1) == 'tga':
                End_coden_Array_right.append(i)
        if len(End_coden_Array_right) < 1 and j < 3:
            End_coden_Array_right.append(int(Seq_len/3 + 20 + End - 1))
        if len(End_coden_Array_right) < 1 and j > 2:
            End_coden_Array_right.append(int(Seq_len/3 + 20 + End - 2))
        begin_two = int(End + 20 + 1)
        end_two = int(Seq_len/3 + Start)
        for m in range(begin_two,end_two):
            end_temp_left1 = TT_seq_Array[m]
            if str(end_temp_left1) == 'tag' or str(end_temp_left1) == 'taa' or str(end_temp_left1) == 'tga':
                End_coden_Array_left.append(m)
        if len(End_coden_Array_left) >= 1:
            begin_two = End_coden_Array_left[len(End_coden_Array_left)-1] + 1
            end_two = int(Seq_len/3 + 20 + End +1)
            for n in range(begin_two,end_two):
                temp_start = TT_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg':
                    Start_coden_Array.append(n)
        if len(End_coden_Array_left) < 1:
            end_three = int(Seq_len/3 + 20 + End + 1)
            for n in range(21 + End,end_three):
                temp_start = TT_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg':
                    Start_coden_Array.append(n)
        if len(Start_coden_Array) < 1 and j > 2:
            Start_coden_Array.append(int(Seq_len/3 + Start - 1))
        if len(Start_coden_Array) < 1 and j < 3:
            Start_coden_Array.append(int(Seq_len/3 + Start))
#################################################################################################################################################################
    if Seq_len % 3 != 0 and cds_end <= Seq_len:
        seq_fold = 7
        begin_one = int(Start + Seq_len)
        end_one = int(Start + Seq_len*2 + 1)
        for i in range(begin_one,end_one):
            end_temp_right1 = ST_seq_Array[i]
            if str(end_temp_right1) == 'tag' or str(end_temp_right1) == 'taa' or str(end_temp_right1) == 'tga':
                End_coden_Array_right.append(i)
        if len(End_coden_Array_right) < 1 and j < 3:
            End_coden_Array_right.append(int(Seq_len + 20 + End - 1))
        if len(End_coden_Array_right) < 1 and j > 2:
            End_coden_Array_right.append(int(Seq_len + 20 + End - 2))
        begin_two = int(End + 20 + 1)
        end_two = int(Seq_len + Start)
        for m in range(begin_two,end_two):
            end_temp_left1 = ST_seq_Array[m]
            if str(end_temp_left1) == 'tag' or str(end_temp_left1) == 'taa' or str(end_temp_left1) == 'tga':
                End_coden_Array_left.append(m)
        if len(End_coden_Array_left) >= 1:
            begin_two = End_coden_Array_left[len(End_coden_Array_left)-1] + 1
            end_two = int(Seq_len + 20 + End + 1)
            for n in range(begin_two,end_two):
                temp_start = ST_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg':
                    Start_coden_Array.append(n)
        if len(End_coden_Array_left) < 1:
            end_three = int(Seq_len + 20 + End + 1)
            for n in range(21 + End,end_three):
                temp_start = ST_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg':
                    Start_coden_Array.append(n)
        if len(Start_coden_Array) < 1 and j > 2:
            Start_coden_Array.append(int(Seq_len + Start - 1))
        if len(Start_coden_Array) < 1 and j < 3:
            Start_coden_Array.append(int(Seq_len + Start))
###################################################################################################################################################################
    if Seq_len % 3 == 0 and cds_end > Seq_len:
        seq_fold = 4
        begin_one = int(Start + Seq_len/3)
        end_one = int(Start + Seq_len/3*2 + 1)
        for i in range(begin_one,end_one):
            end_temp_right1 = FT_seq_Array[i]
            if str(end_temp_right1) == 'tag' or str(end_temp_right1) == 'taa' or str(end_temp_right1) == 'tga':
                End_coden_Array_right.append(i)
        if len(End_coden_Array_right) < 1 and j < 3:
            End_coden_Array_right.append(int(Seq_len/3 + 20 + End - 1))
        if len(End_coden_Array_right) < 1 and j > 2:
            End_coden_Array_right.append(int(Seq_len/3 + 20 + End - 2))
        begin_two = int(End + 20 + 1)
        end_two = int(Seq_len/3 + Start)
        for m in range(begin_two,end_two):
            end_temp_left1 = FT_seq_Array[m]
            if str(end_temp_left1) == 'tag' or str(end_temp_left1) == 'taa' or str(end_temp_left1) == 'tga':
                End_coden_Array_left.append(m)
        if len(End_coden_Array_left) >= 1:
            begin_two = End_coden_Array_left[len(End_coden_Array_left)-1] + 1
            end_two = int(Seq_len/3 + 20 + End + 1)
            for n in range(begin_two,end_two):
                temp_start = FT_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg':
                    Start_coden_Array.append(n)
        if len(End_coden_Array_left) < 1:
            end_three = int(Seq_len/3 + 20 + End + 1)
            for n in range(21 + End,end_three):
                temp_start = FT_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg':
                    Start_coden_Array.append(n)
        if len(Start_coden_Array) < 1 and j > 2:
            Start_coden_Array.append(int(Seq_len/3 + Start - 1))
        if len(Start_coden_Array) < 1 and j < 3:
            Start_coden_Array.append(int(Seq_len/3 + Start))
##########################################################################################################################################################################
    if Seq_len % 3 != 0 and cds_end > Seq_len:
        seq_fold = 8
        begin_one = int(Start + Seq_len)
        end_one = int(Start + Seq_len*2 + 1)
        for i in range(begin_one,end_one):
            end_temp_right1 = ET_seq_Array[i]
            if str(end_temp_right1) == 'tag' or str(end_temp_right1) == 'taa' or str(end_temp_right1) == 'tga':
                End_coden_Array_right.append(i)
        if len(End_coden_Array_right) < 1 and j < 3:
            End_coden_Array_right.append(int(Seq_len + 20 + End - 1))
        if len(End_coden_Array_right) < 1 and j > 2:
            End_coden_Array_right.append(int(Seq_len + 20 + End - 2))
        begin_two = int(End + 20 + 1)
        end_two = int(Seq_len + Start)
        for m in range(begin_two,end_two):
            end_temp_left1 = ET_seq_Array[m]
            if str(end_temp_left1) == 'tag' or str(end_temp_left1) == 'taa' or str(end_temp_left1) == 'tga':
                End_coden_Array_left.append(m)
        if len(End_coden_Array_left) >= 1:
            begin_two = End_coden_Array_left[len(End_coden_Array_left)-1] + 1
            end_two = int(Seq_len + 20 + End + 1)
            for n in range(begin_two,end_two):
                temp_start = ET_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg':
                    Start_coden_Array.append(n)
        if len(End_coden_Array_left) < 1:
            end_three = int(Seq_len + 20 + End + 1)
            for n in range(21 + End,end_three):
                temp_start = ET_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg':
                    Start_coden_Array.append(n)
        if len(Start_coden_Array) < 1 and j > 2:
            Start_coden_Array.append(int(Seq_len + Start - 1))
        if len(Start_coden_Array) < 1 and j < 3:
            Start_coden_Array.append(int(Seq_len + Start))
    return (Start_coden_Array,End_coden_Array_right,seq_fold)
######################################################################
def FindMaxNumber(FrameOrfScore):
    for i in range(len(FrameOrfScore)-1):
        for j in range(i,len(FrameOrfScore)-1):
            if FrameOrfScore[j] > FrameOrfScore[j+1]:
                temp = FrameOrfScore[j]
                FrameOrfScore[j] = FrameOrfScore[j+1]
                FrameOrfScore[j+1] = temp
    return(FrameOrfScore[len(FrameOrfScore)-1])
######################################################################
def Circ_orf_score(Ref_sequence,Quadruple_Ref_sequence,Start_pos_Arr,End_pos,Dir,Index,matrix):
    TempSeuence = InitCodonSeq(0,len(Ref_sequence)-2,3,Ref_sequence)
    TripleSequence = InitCodonSeq(0,len(Quadruple_Ref_sequence)-2,3,Quadruple_Ref_sequence)
    Arr_TempSeuence = TempSeuence.split(' ')
    Arr_TempSeuence.pop()
    Trip_Arr_TempSeuence = TripleSequence.split(' ')
    Trip_Arr_TempSeuence.pop()
    RefCodenSequenceArr = []
    SequenceArr = []
    Frame_orf_Score_arr = []
    Start = ''
    NucleotideSeq_Arr = []
    CodenTotalArr = []
    CodenLengthArr = []
    if int(End_pos) > len(Arr_TempSeuence):
        RefCodenSequenceArr = Trip_Arr_TempSeuence[::]
    else:
        RefCodenSequenceArr = Arr_TempSeuence[::]
    for i in range(0,len(Start_pos_Arr)):
        Cand_Orf_Seq = ''
        frame_orf_score = 0
        temp_start = int(Start_pos_Arr[i])
        temp_stop = int(End_pos)
        for j in range(temp_start,temp_stop):
            Cand_Orf_Seq = Cand_Orf_Seq + RefCodenSequenceArr[j]
        Arr_sequence = SequenceProcessing(Cand_Orf_Seq)
        Arr_coden = InitCodonSeq(0,len(Arr_sequence)-2,3,Arr_sequence)
        Frame_candidate_score_arr = CircSimulation(Arr_coden,len(Arr_sequence)-2,matrix)
        for n in range(len(Frame_candidate_score_arr)):
            frame_orf_score = frame_orf_score + float(Frame_candidate_score_arr[n])
        if frame_orf_score != 0:
            frame_orf_score = frame_orf_score / len(Frame_candidate_score_arr)
        else:
            frame_orf_score = 0
        CodenLengthArr.append(len(Frame_candidate_score_arr))
        Frame_orf_Score_arr.append(frame_orf_score)
        NucleotideSeq_Arr.append(Cand_Orf_Seq)
        CodenTotalArr.append(Arr_coden)
    frame_max_score = FindMaxNumber(Frame_orf_Score_arr)
    for w in range(len(Frame_orf_Score_arr)):
        if float(frame_max_score) == Frame_orf_Score_arr[w]:
            Start = Start_pos_Arr[w]
    return (frame_max_score,Start,NucleotideSeq_Arr[w],CodenTotalArr[w],CodenLengthArr[w])
#########################################
def Small_orf_score_max(TT_seq,FT_seq,ST_seq,ET_seq,Start_pos_Arr,End_pos_Arr,matrix,six,End_mlcds,Seq_len,score,Start_mlcds,sequence,length):
    cds_end = End_mlcds * 3 + 20 * 3 + six
    if Seq_len % 3 == 0 and cds_end <= Seq_len:
        Ref_sequence = TT_seq
    if Seq_len % 3 != 0 and cds_end <= Seq_len:
        Ref_sequence = ST_seq
    if Seq_len % 3 == 0 and cds_end > Seq_len:
        Ref_sequence = FT_seq
    if Seq_len % 3 != 0 and cds_end > Seq_len:
        Ref_sequence = ET_seq
    Arr_TempSeuence = Ref_sequence.split(' ')
    Arr_TempSeuence.pop()
    RefCodenSequenceArr = []
    SequenceArr = []
    frame_score = []
    Start = []
    End = []
    NucleotideSeq_Arr = []
    CodenLengthArr = []
    Starts = []
    pp = 0
    if Start_pos_Arr[0] < End_pos_Arr[len(End_pos_Arr)-1]:
        for i in range(len(End_pos_Arr)):
            for j in range(len(Start_pos_Arr)):
                if Start_pos_Arr[j] > pp and Start_pos_Arr[j] < End_pos_Arr[i]:
                    Cand_Orf_Seq = ''
                    frame_orf_score = 0
                    temp_start = int(Start_pos_Arr[j])
                    temp_stop = int(End_pos_Arr[i])
                    Start1 = temp_start*3 + six +1
                    End1 = temp_stop*3 + six + 3
                    Start.append(int(Start1))
                    End.append(int(End1))
                    for m in range(temp_start,temp_stop+1):
                        Cand_Orf_Seq = Cand_Orf_Seq + Arr_TempSeuence[m]
                    Arr_sequence = SequenceProcessing(Cand_Orf_Seq)
                    Arr_coden = InitCodonSeq(0,len(Arr_sequence)-2,3,Arr_sequence)
                    frame_orf_score = CircSimulationmlcds(Arr_coden,matrix)
                    CodenLengthArr.append(len(Arr_sequence))
                    frame_orf_score = frame_orf_score[0]/len(Arr_sequence)
                    frame_score.append(frame_orf_score)
                    NucleotideSeq_Arr.append(Cand_Orf_Seq)
            pp = End_pos_Arr[i]
    else:
        frame_score.append(score)
        Start.append(Start_mlcds)
        End.append(End_mlcds)
        NucleotideSeq_Arr.append(sequence)
        CodenLengthArr.append(length)
    return (frame_score,Start,End,NucleotideSeq_Arr,CodenLengthArr)
##########################################
def CircSimulationmlcds(line_sequence, hash_score):
    mlcds_score = []
    Sub_Frame_orf_seq = line_sequence.split(' ')
    Sub_Frame_orf_seq.pop()
    number = 0
    for n in range(0, len(Sub_Frame_orf_seq) - 1):
        temp1 = Sub_Frame_orf_seq[n] + Sub_Frame_orf_seq[n + 1]
        num_temp = re.compile('[atcg]{6}')
        if num_temp.match(temp1):
            number = float(number) + float(hash_score[temp1])
    mlcds_score.append(number)
    return (mlcds_score)
##########################################
def mainProcess(inputFile, codonArr, hash_matrix):
    label_Arr, FastA_seq_Arr = Reading(inputFile)
    for i in range(len(label_Arr)):
        Label = label_Arr[i]
        Seq = FastA_seq_Arr[i]
        re_Seq = Seq[::-1]
        Double_Seq = Seq + Seq
        Re_Double_Seq = re_Seq + re_Seq
        Three_seq = Seq + Seq + Seq
        Four_seq = Seq + Seq + Seq + Seq
        Re_Three_seq = re_Seq + re_Seq + re_Seq
        Re_Four_seq = re_Seq + re_Seq + re_Seq + re_Seq
        Seven_seq = Seq + Seq + Seq + Seq + Seq + Seq + Seq
        Eight_seq = Seq + Seq + Seq + Seq + Seq + Seq + Seq + Seq
        Re_Seven_seq = re_Seq + re_Seq + re_Seq + re_Seq + re_Seq + re_Seq + re_Seq
        Re_Eight_seq = re_Seq + re_Seq + re_Seq + re_Seq + re_Seq + re_Seq + re_Seq + re_Seq
        Seq_Arr = SequenceProcessing(Seq)
        cis_seq_Arr = SequenceProcessing(Double_Seq)
        Re_dou_Arr = SequenceProcessing(Re_Double_Seq)
        Three_ref_arr = SequenceProcessing(Three_seq)
        Re_Three_ref_arr = SequenceProcessing(Re_Three_seq)
        Four_ref_arr = SequenceProcessing(Four_seq)
        Re_Four_ref_arr = SequenceProcessing(Re_Four_seq)
        Seven_ref_arr = SequenceProcessing(Seven_seq)
        Re_Seven_ref_arr = SequenceProcessing(Re_Seven_seq)
        Eight_ref_arr = SequenceProcessing(Eight_seq)
        Re_Eight_ref_arr = SequenceProcessing(Re_Eight_seq)
        RnaOrfSequence = []
        RnaOrfStart = []
        RnaOrfStop = []
        RnaOrfLength = []
        RnaOrfScore = []
        RnaOrfFrameIndex = []
        RnaOrfDirectory = []
        CDS_Score = []
        CDS_length = []
        CDS_conservation = []
        CDS_sequence = []
        ML_RL = []
        Sub_Frame = []
        Seq_len = len(cis_seq_Arr)
        if Seq_len/2 > 61:
            for j in range(0, 6):
                TempStr = ''
                if j < 3 :
                    TempStr = InitCodonSeq(j,Seq_len-2,3,cis_seq_Arr)
                if 2 < j < 6 :
                    TempStr = InitCodonSeq(j-3,Seq_len-2,3,Re_dou_Arr)
                Each_frame_score = CircSimulation(TempStr,hash_matrix)
                Sub_Frame_orf_Score,Sub_Frame_orf,Sub_Start,Sub_End = MaxSubseqPlus(Each_frame_score,TempStr)
                Sub_Frame_orf_seq = Sub_Frame_orf.split(' ')
                Sub_Frame_orf_seq.pop()
                MLCDS_length = len(Sub_Frame_orf_seq) * 3
                CDS_seq = "".join(Sub_Frame_orf_seq)
                MLSCDS_score = CircSimulationmlcds(Sub_Frame_orf,hash_matrix)
                MLSCDS_score_mean = MLSCDS_score[0] / MLCDS_length
                Sub_Frame.append([MLSCDS_score_mean,CDS_seq,MLCDS_length,Sub_Start,Sub_End,j])
            MLSCDS_score_mean = list(map(lambda x:x[0],Sub_Frame))
            max_index = MLSCDS_score_mean.index(max(MLSCDS_score_mean))
            MLSCDS_score_mean,CDS_seq,MLCDS_length,Sub_Start,Sub_End,j = Sub_Frame[max_index]
            Coding_method_j = j
            if j <= 2:
                score, number = add_m6a_info(Seq)
                feature_score = get_feature_score(Seq)
            else:
                score, number = add_m6a_info(re_Seq)
                feature_score = get_feature_score(re_Seq)
            lengthbi = MLCDS_length / Seq_len
            Frame_start_coden_sequence = []
            Frame_stop_coden_sequence = []
            Directory = ''
            if j < 3:
                Directory = '+'
            else:
                Directory = '-'
            CDS_Score.append(MLSCDS_score_mean)
            CDS_length.append(MLCDS_length)
            CDS_sequence.append(CDS_seq)
            ML_RL.append(lengthbi)
            Seq_len = len(Seq_Arr)
            if j < 3 :
                Three_ref_arr = InitCodonSeq(j,Seq_len*3-2,3,Three_ref_arr)
                Four_ref_arr = InitCodonSeq(j,Seq_len*4-2,3,Four_ref_arr)
                Seven_ref_arr = InitCodonSeq(j,Seq_len*7-2,3,Seven_ref_arr)
                Eight_ref_arr = InitCodonSeq(j,Seq_len*8-2,3,Eight_ref_arr)
            if 2 < j < 6 :
                Three_ref_arr = InitCodonSeq(j,Seq_len*3-2,3,Re_Three_ref_arr)
                Four_ref_arr = InitCodonSeq(j,Seq_len*4-2,3,Re_Four_ref_arr)
                Seven_ref_arr = InitCodonSeq(j,Seq_len*7-2,3,Re_Seven_ref_arr)
                Eight_ref_arr = InitCodonSeq(j,Seq_len*8-2,3,Re_Eight_ref_arr)
            Frame_start_coden_sequence, Frame_stop_coden_sequence, seq_fold= FindStartCodne(Three_ref_arr,Four_ref_arr,Seven_ref_arr,Eight_ref_arr,Sub_Start,Sub_End,j,Seq_len)
            if len(Frame_start_coden_sequence) >= 1 and len(Frame_stop_coden_sequence) >= 1:
                RnaOrfScore, RnaOrfStart, RnaOrfStop, RnaOrfSequence, RnaOrfLength = Small_orf_score_max(Three_ref_arr,Four_ref_arr,Seven_ref_arr,Eight_ref_arr,Frame_start_coden_sequence,Frame_stop_coden_sequence,hash_matrix,j,Sub_End,Seq_len,MLSCDS_score_mean,Sub_Start,CDS_seq,MLCDS_length)
                RnaOrfFrameIndex.append(j)
                RnaOrfDirectory.append(Directory)
            Str_UP = CDS_seq.upper()
            H_temp = '>' + 'Human' + '|'
            H_Cov_str_temp = H_temp + str(Label) + '|' + str(j) + '\n'
            H_Cov_str = H_Cov_str_temp + Str_UP + '\n'
            CDS_file = Circ_Dir + '/' + str(Label) + '|' + str(j)
            CDS_path_name_temp = '\t' + CDS_file + '\t'
            CDS_path_name = CDS_path_name_temp.replace('\t', "'")
            CDS_path = open(CDS_file, 'w')
            CDS_path.write(H_Cov_str)
            CDS_path.close()
            Conservation_Score_Read = os.popen(
                'PhyloCSF 100vertebrates --strategy=fixed --minCodons=30 --frames=3 --removeRefGaps  ' + CDS_path_name + ' |awk \'{print $3}\'').read()
            Conservation_Score_Read = Conservation_Score_Read.replace('\n', '')
            CDS_conservation.append(Conservation_Score_Read)
            for m in range(len(RnaOrfStart)):
                RnaOutInformation = ''
                RnaOutInformation = RnaOutInformation + 'Index:' + str(RnaOrfFrameIndex[0]) + '\t'
                RnaOutInformation = RnaOutInformation + 'Directory:' + str(RnaOrfDirectory[0]) + '\t'
                RnaOutInformation = RnaOutInformation + 'Sequence:' + str(RnaOrfSequence[m]) + '\t'
                RnaOutInformation = RnaOutInformation + 'Sequence_fold:' + str(seq_fold) + '\t'
                RnaOutInformation = RnaOutInformation + 'Score:' + str(RnaOrfScore[m]) + '\t'
                RnaOutInformation = RnaOutInformation + 'Start:' + str(RnaOrfStart[m]) + '\t'
                RnaOutInformation = RnaOutInformation + 'End:' + str(RnaOrfStop[m]) + '\t'
                RnaOutInformation = RnaOutInformation + 'CDS_Score:' + str(MLSCDS_score_mean) + '\t'
                RnaOutInformation = RnaOutInformation + 'Length:' + str(RnaOrfLength[m]) + '\t'
                RnaOutInformation = RnaOutInformation + 'origin_Seq_len:' + str(Seq_len) + '\t'
                RnaOutInformation = RnaOutInformation + 'Coding_method_j:' + str(Coding_method_j) + '\t'
                RnaOutInformation = RnaOutInformation + 'circRNA sequence:' + str(Seq) + '\t'
                RnaOutInformation = 'CircRNA_ID:' + str(Label) + '\t' + RnaOutInformation + '\n'
                Temp_Out_Information1.write(RnaOutInformation)
            CDS_Feature = ''
            CDS_Feature = CDS_Feature + str(CDS_conservation[0]) + ' '
            CDS_Feature = CDS_Feature + str(CDS_Score[0]) + ' '
            CDS_Feature = CDS_Feature + str(CDS_length[0]) + ' '
            CDS_Feature = CDS_Feature + str(ML_RL[0]) + ' '
            CDS_Feature = CDS_Feature + str(number) + ' '
            CDS_Feature = CDS_Feature + str(feature_score) + ' '
            CDS_Feature = CDS_Feature + str(CDS_sequence[0])
            CDS_Feature = str(Label) + '|' + str(j) + ' ' + CDS_Feature + '\n'
            Temp_Out_Feature1.write(CDS_Feature)
            if i == len(label_Arr)-1:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}")
            else:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}\r\n")
        elif 3 <= Seq_len/2 <= 61:
            Seq_len = int(Seq_len/2)
            RnaOrfSequence.append(FastA_seq_Arr[i])
            RnaOrfStart.append('1')
            RnaOrfStop.append(Seq_len)
            RnaOrfLength.append(Seq_len)
            RnaOrfFrameIndex.append('null')
            RnaOrfDirectory.append('null')
            TempStr = InitCodonSeq(0, Seq_len - 2, 3, Seq_Arr)
            score, number = add_m6a_info(Seq)
            feature_score = get_feature_score(Seq)
            FrameScore = CircSimulationmlcds(TempStr, hash_matrix)
            MLSCDS_score_mean = FrameScore[0]/Seq_len
            TempArray = TempStr.split(' ')
            TempArray.pop()
            CDS_seq = "".join(TempArray)
            RnaOrfScore.append(MLSCDS_score_mean)
            CDS_Score.append(MLSCDS_score_mean)
            CDS_length.append(Seq_len)
            CDS_conservation.append('0')
            CDS_sequence.append(FastA_seq_Arr[i])
            RnaOutInformation = ''
            RnaOutInformation = RnaOutInformation + 'Index:' + str(RnaOrfFrameIndex[0])  + '\t'
            RnaOutInformation = RnaOutInformation + 'Directory:' + str(RnaOrfDirectory[0]) + '\t'
            RnaOutInformation = RnaOutInformation + 'Sequence:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'Sequence_fold:' + '0' + '\t'
            RnaOutInformation = RnaOutInformation + 'Score:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'Start:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'End:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'CDS_Score:' + str(MLSCDS_score_mean) + '\t'
            RnaOutInformation = RnaOutInformation + 'Length:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'origin_Seq_len:' + str(Seq_len) + '\t'
            RnaOutInformation = RnaOutInformation + 'Coding_method_j:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'circRNA sequence:' + str(Seq) + '\t'
            RnaOutInformation = 'circRNA_ID:' + str(Label) + '\t' +  RnaOutInformation + '\n'
            Temp_Out_Information2.write(RnaOutInformation)
            CDS_Feature = ''
            CDS_Feature = CDS_Feature + str(CDS_conservation[0]) + ' '
            CDS_Feature = CDS_Feature + 'null' + ' '
            CDS_Feature = CDS_Feature + 'null' + ' '
            CDS_Feature = CDS_Feature + str('0') + ' '
            CDS_Feature = CDS_Feature + str(number) + ' '
            CDS_Feature = CDS_Feature + str(feature_score) + ' '
            CDS_Feature = CDS_Feature + 'null'
            CDS_Feature = str(Label) + '|' + str('0') + ' ' + CDS_Feature + '\n'
            Temp_Out_Feature2.write(CDS_Feature)
            if i == len(label_Arr) - 1:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}")
            else:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}\r\n")
        else:
            Seq_len = int(Seq_len/2)
            n = 0
            RnaOutInformation = ''
            RnaOutInformation = RnaOutInformation + 'Index:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'Directory:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'Sequence:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'Sequence_fold:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'Score:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'Start:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'End:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'CDS_Score:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'Length:' + 'null' + '\t'
            RnaOutInformation = RnaOutInformation + 'origin_Seq_len:' + str(Seq_len) + '\t'
            RnaOutInformation = RnaOutInformation + 'Coding_method_j:' + 'null'+ '\t'
            RnaOutInformation = RnaOutInformation + 'circRNA sequence:' + str(Seq) + '\t'
            RnaOutInformation = 'circRNA_ID:' + str(Label) + '\t' +  RnaOutInformation + '\n'
            Temp_Out_Information3.write(RnaOutInformation)
            CDS_Feature = ''
            CDS_Feature = CDS_Feature + 'null' + ' '
            CDS_Feature = CDS_Feature + 'null' + ' '
            CDS_Feature = CDS_Feature + 'null' + ' '
            CDS_Feature = CDS_Feature + 'null' + ' '
            CDS_Feature = CDS_Feature + 'null' + ' '
            CDS_Feature = CDS_Feature + 'null' + ' '
            CDS_Feature = CDS_Feature + 'null'
            CDS_Feature = str(Label) + '|' + str('0') + ' ' + CDS_Feature + '\n'
            Temp_Out_Feature3.write(CDS_Feature)
            if i == len(label_Arr) - 1:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}")
            else:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}\r\n")
######################################################################
def add_m6a_info(seq):
    with open('{}/seq_sramp.fa'.format(Circ_Dir),'w+') as seq_sramp:
        seq = '>seq\n' + seq
        seq_sramp.write(seq)
    uuid_fa = str(uuid.uuid1()) + "_temp_m6a.txt"
    os.system("perl runsramp.pl %s/seq_sramp.fa  %s full" % (Circ_Dir,uuid_fa))
    if not os.path.isfile(uuid_fa):
        return [0,0]
    os.rename(uuid_fa,os.path.join(Circ_Dir,'temp_m6a.txt'))
    input = pd.read_csv("{}/temp_m6a.txt".format(Circ_Dir), sep="\t")[['Seq_ID', 'Position', 'Score(Combined)', 'Classification']]
    input_m6a = input[input.Classification.apply(lambda r: "Non" not in r)]
    if input_m6a.empty: 
        return [0,0]
    else:
        input_m6a_group = input_m6a.groupby(["Seq_ID", "Classification"]).agg({
            'Position': len,
            'Score(Combined)': sum
        }).reset_index()
        input_m6a_group['Score(Combined)'] = input_m6a_group.loc[input_m6a_group.Classification.str.contains("Low|Very high|High|Moderate"),'Score(Combined)']
        input_m6a_group = input_m6a_group[['Seq_ID', 'Position', 'Score(Combined)']]
        input_m6a_group.columns = ['Seq_ID', 'Count', 'Score']
        m6a_group_final = input_m6a_group.groupby(["Seq_ID"]).sum()
        m6a_group_final = m6a_group_final[['Count', 'Score']].values / len(seq)
        return m6a_group_final[0]
##############################################################################
def get_feature_score(seq):
    with open('{}/seq_featrue_score.fa'.format(Circ_Dir),'w+') as seq_featrue_score:
        seq = '>seq\n' + seq
        seq_featrue_score.write(seq)
    importr("LncFinder")
    r('''
    getFeatureScore <- function(){{
        Seqs <- seqinr::read.fasta("{}")
        score <- compute_FickettScore(Seqs, label = NULL, on.ORF = TRUE, auto.full = TRUE, parallel.cores = 2)
        score$Fickett.Score
    }}
    '''.format(os.path.join(os.getcwd(),Circ_Dir,'seq_featrue_score.fa')))
    feature_score = r['getFeatureScore']()
    return list(feature_score)[0]

if __name__ == '__main__':
    usage = "[options]"
    parse = OptionParser(usage)
    parse.add_option('-f', '--file', dest='file', action='store', metavar='input files',
                     help='input files' 
                     '(Required) Input transcript nucleotide sequences')
    parse.add_option('-o', '--out', dest='outfile', action='store', metavar='output files path',
                     help='output files path'
                     '(Required) The path where you want to output')
    parse.add_option('-t', '--top', dest='top', action='store', metavar='top rank protein products', default='all',
                     help='top rank protein products'
                     '(Optional) The top rank protein products you want')
    (options, args) = parse.parse_args()
    inPutFileName = options.file
    outPutFileName = options.outfile
    top_num = options.top
    curPath = cur_file_dir()
    MatrixPath = curPath + "/CICADA_Parameters/circRNA_Score_Matrix"
    inMatrix = open(MatrixPath)
    Matrix = inMatrix.read()
    inMatrix.close()
    Circ_Dir = outPutFileName + '_Tmp_Dir'
    subprocess.call('mkdir ' + Circ_Dir + '', shell=True)
    Out_Information1 = Circ_Dir + '/Information1'
    Out_Information2 = Circ_Dir + '/Information2'
    Out_Information3 = Circ_Dir + '/Information3'
    Out_Information = Circ_Dir + '/Information'
    Temp_Out_Information1 = open(Out_Information1, 'w')
    Temp_Out_Information2 = open(Out_Information2, 'w')
    Temp_Out_Information3 = open(Out_Information3, 'w')
    Out_Feature1 = Circ_Dir + '/Feature1'
    Out_Feature2 = Circ_Dir + '/Feature2'
    Out_Feature3 = Circ_Dir + '/Feature3'
    Out_Feature = Circ_Dir + '/Feature'
    Temp_Out_Feature1 = open(Out_Feature1, 'w')
    Temp_Out_Feature2 = open(Out_Feature2, 'w')
    Temp_Out_Feature3 = open(Out_Feature3, 'w')
    Out_Feature_fasta = Out_Feature + "_fasta.fa"
    Temp_Out_Feature_fasta = open(Out_Feature_fasta, 'w')
    CDS_path = Circ_Dir + '/C_score_write'
    Alphabet = ['ttt', 'ttc', 'tta', 'ttg', 'tct', 'tcc', 'tca', 'tcg', 'tat', 'tac', 'tgt', 'tgc', 'tgg', 'ctt', 'ctc',
                'cta', 'ctg', 'cct', 'ccc', 'cca', 'ccg', 'cat', 'cac', 'caa', 'cag', 'cgt', 'cgc', 'cga', 'cgg', 'att',
                'atc', 'ata', 'atg', 'act', 'acc', 'aca', 'acg', 'aat', 'aac', 'aaa', 'aag', 'agt', 'agc', 'aga', 'agg',
                'gtt', 'gtc', 'gta', 'gtg', 'gct', 'gcc', 'gca', 'gcg', 'gat', 'gac', 'gaa', 'gag', 'ggt', 'ggc', 'gga',
                'ggg']
    Matrix_hash = {}
    Matrix_Arr = Matrix.split('\n')
    length = len(Matrix_Arr) - 1
    del Matrix_Arr[length]
    for line in Matrix_Arr:
        each = line.split('\t')
        key = each[0]
        value = each[1]
        Matrix_hash[key] = value
    inFiles = open(inPutFileName)
    inFilesArr = inFiles.read()
    inFileNum = inFilesArr.split('\n')
    inFileLen = len(inFileNum) - 1
    inFiles.close()
    sequence_Arr = inFilesArr.split('\n')
    ARRAY = TwoLineFasta(sequence_Arr)
    Label_Array, FastA_Seq_Array = Tran_checkSeq(ARRAY)
    inFileLength = len(Label_Array)
    TOT_STRING = []
    for i in range(len(Label_Array)):
        tmp_label_one = Label_Array[i]
        tmp_label = tmp_label_one.replace('\r', '')
        tmp_seq = FastA_Seq_Array[i]
        Temp_Seq = tmp_seq.replace('\r', '')
        TOT_STRING.append(tmp_label)
        TOT_STRING.append(Temp_Seq)
    mainProcess(TOT_STRING, Alphabet, Matrix_hash)
    Temp_Out_Information1.close()
    Temp_Out_Information2.close()
    Temp_Out_Information3.close()
    Temp_Out_Feature1.close()
    Temp_Out_Feature2.close()
    Temp_Out_Feature3.close()
    Temp_Out_Feature_fasta.close()
    Pre_mian(Out_Feature1,Out_Feature2,Out_Feature3, Out_Information1, Out_Information2, Out_Information3, top_num)
    os.remove(Circ_Dir + '/Feature1')
    os.remove(Circ_Dir + '/Feature2')
    os.remove(Circ_Dir + '/Feature3')
    os.remove(Circ_Dir + '/Information1')
    os.remove(Circ_Dir + '/Information2')
    os.remove(Circ_Dir + '/Information3')