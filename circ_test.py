###################################################################### -*-coding : utf-8-*-
import re
import os
import joblib
import pandas as pd
import optparse
import time
import subprocess
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import sys
import faulthandler;faulthandler.enable()
# run shell
#export R_HOME=/home/software/R-3.5.0
#export LD_LIBRARY_PATH=/home/software/R-3.5.0/lib

R_HOME = '/home/ZXY/software/R-3.5.0'
LD_LIBRARY_PATH = '/home/ZXY/software/R-3.5.0/lib'
os.environ['R_HOME'] = R_HOME
#
#
#if LD_LIBRARY_PATH not in os.environ['LD_LIBRARY_PATH']:
#    os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH
#    os.execv('/usr/bin/python3', [","] + sys.argv)
##########################################################################################
def xgboostPredict(data):
    xgb = joblib.load('/home/ZXY/peptide/circRNA/model_own/hh_cir_xgb_test_Predict.model')
    pre = xgb.predict_proba(data)
    return pre


def rfPredict(data):
    rf = joblib.load('/home/ZXY/peptide/circRNA/jiaocha/model_saved/Cir_RF_Predict.model')
    pre = rf.predict_proba(data)
    return pre


def Pre_mian(fname,Out_Information, model):
    start = time.time()
    data = pd.read_csv(fname,sep=' ',header=None)
    dataused = data.iloc[:,[1,2,3,4,5,7]].values
    # list: [1. 1. 1. ... 1. 1. 0.]
    if model == 'xgb':
        pre = xgboostPredict(dataused)
    elif model == 'rf':
        pre = rfPredict(dataused)
    else:
        print("Please input correct model name.")
    data[10] = pd.Series(pre[:,0])
    data[11] = pd.Series(pre[:,1])
    data.to_csv(fname,sep='\t', index=False, header=False)
    id_dict = dict(map(lambda x:[x[0].split('|')[0],(x[10],x[11])],data.values.tolist()))
    information = pd.read_csv(Out_Information,sep='\t',header=None)
    information[9] = information[0].apply(lambda x:id_dict.get(x.split(':')[1])[0])
    information[10] = information[0].apply(lambda x:id_dict.get(x.split(':')[1])[1])
    information.to_csv(Out_Information,sep='\t',index=False,header=False)
    end = time.time()
#over######################################################################
def cur_file_dir():
    path = sys.path[0]
    if os.path.isdir(path):
       return path
    elif os.path.isfile(path):
       return os.path.dirname(path)
#over######################################################################
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
#over######################################################################
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
    #LogResult = Temp_Log + '.log'
    #LOG_FILE = open(LogResult,'w')
    num = 0
    for i in range(len(label_Arr)):
        Label = label_Arr[num]
        Seq = FastA_seq_Arr[num]
        tran_fir_seq = Seq.lower()
        tran_sec_seq_one = tran_fir_seq.replace('u','t')
        tran_sec_seq = tran_sec_seq_one.replace('\r','')
        if 'n' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (n),please checkout your sequence again' + '\n'
            #LOG_FILE.write(LogString)
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'w' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (w),please checkout your sequence again' + '\n'
            #LOG_FILE.write(LogString)
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'd' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (d),please checkout your sequence again' + '\n'
            #LOG_FILE.write(LogString)
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'r' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (r),please checkout your sequence again' + '\n'
            #LOG_FILE.write(LogString)
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 's' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (s),please checkout your sequence again' + '\n'
            #LOG_FILE.write(LogString)
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'y' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (y),please checkout your sequence again' + '\n'
            #LOG_FILE.write(LogString)
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        if 'm' in tran_sec_seq:
            LogString = Label + ' ' + 'contain unknow nucleotide (m),please checkout your sequence again' + '\n'
            #LOG_FILE.write(LogString)
            del label_Arr[num]
            del FastA_seq_Arr[num]
            continue
        num = int(num) + int(1)
    #LOG_FILE.close()
    return (label_Arr,FastA_seq_Arr)
#over######################################################################
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
#over######################################################################
def SequenceProcessing(sequence):
    length = len(sequence)
    tran_lower_seq = sequence.lower()
    tran_cis_seq = tran_lower_seq.replace('u','t')
    cis_sequence_Arr = list(tran_cis_seq)
    return cis_sequence_Arr
#over######################################################################
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
#over######################################################################
def CircSimulation(line_sequence,hash_score):
    CodonScore = []# coden score for every frame
    #Total_str_ARR.append(TempStr)
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
                number = float(number) + float(hash_score[temp1])### sorce matrix is very imporant
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
######################################################################3倍—7倍-4倍-8倍
def FindStartCodne(TT_seq,FT_seq,ST_seq,ET_seq,Start,End,j,Seq_len):
    TT_seq_Array = TT_seq.split(' ')###coden array
    TT_seq_Array.pop()
    ST_seq_Array = ST_seq.split(' ')###coden array
    ST_seq_Array.pop()
    FT_seq_Array = FT_seq.split(' ')###coden array
    FT_seq_Array.pop()
    ET_seq_Array = ET_seq.split(' ')###coden array
    ET_seq_Array.pop()
    Start_coden_Array = []
    End_coden_Array_left = []
    End_coden_Array_right = []
    cds_end = End * 3 + 20 * 3 + j
    if Seq_len % 3 == 0 and cds_end <= Seq_len:
        begin_one = int(Start + Seq_len/3)
        end_one = int(Seq_len/3 * 2 + Start + 1)
        for i in range(begin_one,end_one):
            end_temp_right1 = TT_seq_Array[i]
            if str(end_temp_right1) == 'tag' or str(end_temp_right1) == 'taa' or str(end_temp_right1) == 'tga':
                End_coden_Array_right.append(i)
        begin_two = int(End + 20 + 1)
        end_two = int(Seq_len/3 + Start)
        for m in range(begin_two,end_two):
            end_temp_left1 = TT_seq_Array[m]
            if str(end_temp_left1) == 'tag' or str(end_temp_left1) == 'taa' or str(end_temp_left1) == 'tga':
                End_coden_Array_left.append(m)
        if len(End_coden_Array_left) >= 1:
            End_coden_left = End_coden_Array_left[len(End_coden_Array_left)-1] + 1
            end_three = int(Seq_len/3 + 20 + End +1)
            for n in range(End_coden_left,end_three):
                temp_start = TT_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg': # possibility can find other feature in sequence in fornt of start coden
                    Start_coden_Array.append(n)
        else:
            end_three = int(Seq_len/3 + 20 + End + 1)
            for n in range(21 + End,end_three):
                temp_start = TT_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg': # possibility can find other feature in sequence in fornt of start coden
                    Start_coden_Array.append(n)
#################################################################################################################################################################
    if Seq_len % 3 != 0 and cds_end <= Seq_len:
        begin_one = int(Start + Seq_len)
        end_one = int(Start + Seq_len*2 + 1)
        for i in range(begin_one,end_one):
            end_temp_right1 = ST_seq_Array[i]
            if str(end_temp_right1) == 'tag' or str(end_temp_right1) == 'taa' or str(end_temp_right1) == 'tga':
                End_coden_Array_right.append(i)
        begin_two = int(End + 20 + 1)
        end_two = int(Seq_len + Start)
        for m in range(begin_two,end_two):
            end_temp_left1 = ST_seq_Array[m]
            if str(end_temp_left1) == 'tag' or str(end_temp_left1) == 'taa' or str(end_temp_left1) == 'tga':
                End_coden_Array_left.append(m)
        if len(End_coden_Array_left) >= 1:
            End_coden_left = End_coden_Array_left[len(End_coden_Array_left)-1] + 1
            end_three = int(Seq_len + 20 + End + 1)
            for n in range(End_coden_left,end_three):
                temp_start = ST_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg': # possibility can find other feature in sequence in fornt of start coden
                    Start_coden_Array.append(n)
        else:
            end_three = int(Seq_len + 20 + End + 1)
            for n in range(21 + End,end_three):
                temp_start = ST_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg': # possibility can find other feature in sequence in fornt of start coden
                    Start_coden_Array.append(n)
###################################################################################################################################################################
    if Seq_len % 3 == 0 and cds_end > Seq_len:
        begin_one = int(Start + Seq_len/3)
        end_one = int(Start + Seq_len/3*2 + 1)
        for i in range(begin_one,end_one):
            end_temp_right1 = FT_seq_Array[i]
            if str(end_temp_right1) == 'tag' or str(end_temp_right1) == 'taa' or str(end_temp_right1) == 'tga':
                End_coden_Array_right.append(i)
        begin_two = int(End + 20 + 1)
        end_two = int(Seq_len/3 + Start)
        for m in range(begin_two,end_two):
            end_temp_left1 = FT_seq_Array[m]
            if str(end_temp_left1) == 'tag' or str(end_temp_left1) == 'taa' or str(end_temp_left1) == 'tga':
                End_coden_Array_left.append(m)
        if len(End_coden_Array_left) >= 1:
            End_coden_left = End_coden_Array_left[len(End_coden_Array_left)-1] + 1
            end_three = int(Seq_len/3 + 20 + End + 1)
            for n in range(End_coden_left,end_three):
                temp_start = FT_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg': # possibility can find other feature in sequence in fornt of start coden
                    Start_coden_Array.append(n)
        else:
            end_three = int(Seq_len/3 + 20 + End + 1)
            for n in range(21 + End,end_three):
                temp_start = FT_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg': # possibility can find other feature in sequence in fornt of start coden
                    Start_coden_Array.append(n)
##########################################################################################################################################################################
    if Seq_len % 3 != 0 and cds_end > Seq_len:
        begin_one = int(Start + Seq_len)
        end_one = int(Start + Seq_len*2 + 1)
        for i in range(begin_one,end_one):
            end_temp_right1 = ET_seq_Array[i]
            if str(end_temp_right1) == 'tag' or str(end_temp_right1) == 'taa' or str(end_temp_right1) == 'tga':
                End_coden_Array_right.append(i)
        begin_two = int(End + 20 + 1)
        end_two = int(Seq_len + Start)
        for m in range(begin_two,end_two):
            end_temp_left1 = ET_seq_Array[m]
            if str(end_temp_left1) == 'tag' or str(end_temp_left1) == 'taa' or str(end_temp_left1) == 'tga':
                End_coden_Array_left.append(m)
        if len(End_coden_Array_left) >= 1:
            End_coden_left = End_coden_Array_left[len(End_coden_Array_left)-1] + 1
            end_three = int(Seq_len + 20 + End + 1)
            for n in range(End_coden_left,end_three):
                temp_start = ET_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg': # possibility can find other feature in sequence in fornt of start coden
                    Start_coden_Array.append(n)
        else:
            end_three = int(Seq_len + 20 + End + 1)
            for n in range(21 + End,end_three):
                temp_start = ET_seq_Array[n]
                if str(temp_start) == 'atg' or str(temp_start) == 'ttg' or str(temp_start) == 'ctg' or str(temp_start) == 'gtg': # possibility can find other feature in sequence in fornt of start coden
                    Start_coden_Array.append(n)
    return (Start_coden_Array,End_coden_Array_right)
#def PhyloCSF(Write_Path,Read_Path):
    #os.system('/linuxdata3/zxy/software/PhyloCSF-master/PhyloCSF 100vertebrates --files '+Write_Path+' --strategy=fixed --minCodons=30 --removeRefGaps > '+Read_Path+'')
######################################################################
def FindMaxNumber(FrameOrfScore):
    for i in range(len(FrameOrfScore)-1):
        for j in range(i,len(FrameOrfScore)-1):
            if FrameOrfScore[j] > FrameOrfScore[j+1]:
                temp = FrameOrfScore[j]
                FrameOrfScore[j] = FrameOrfScore[j+1]
                FrameOrfScore[j+1] = temp
    return(FrameOrfScore[len(FrameOrfScore)-1])
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
def Small_orf_score_max(TT_seq,FT_seq,ST_seq,ET_seq,Start_pos_Arr,End_pos_Arr,matrix,j,End,Seq_len):
    cds_end = End * 3 + 20 * 3 + j
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
    Frame_orf_Score_arr = []
    Start = ''
    End = ''
    NucleotideSeq_Arr = []
    CodenLengthArr = []
    Starts = []
    pp = 0
    if Start_pos_Arr[len(Start_pos_Arr)-1] < End_pos_Arr[0]:
        for i in range(len(End_pos_Arr)):
            for j in range(len(Start_pos_Arr)):
                if Start_pos_Arr[j] > pp and Start_pos_Arr[j] < End_pos_Arr[i]:
                    Cand_Orf_Seq = ''
                    frame_orf_score = 0
                    temp_start = int(Start_pos_Arr[j])
                    temp_stop = int(End_pos_Arr[i])
                    for m in range(temp_start,temp_stop+1):
                        Cand_Orf_Seq = Cand_Orf_Seq + Arr_TempSeuence[m]
                    Arr_sequence = SequenceProcessing(Cand_Orf_Seq)
                    Arr_coden = InitCodonSeq(0,len(Arr_sequence)-2,3,Arr_sequence)
                    frame_orf_score = CircSimulationmlcds(Arr_coden,matrix)
                    CodenLengthArr.append(len(Arr_sequence))
                    frame_orf_score = frame_orf_score[0]/len(Arr_sequence)
                    Frame_orf_Score_arr.append(frame_orf_score)
                    NucleotideSeq_Arr.append(Cand_Orf_Seq)
                    Starts.append(Start_pos_Arr[j])
            pp = End_pos_Arr[i]
        frame_max_score = FindMaxNumber(Frame_orf_Score_arr)
        for w in range(len(Starts)):
            if float(frame_max_score) == float(Frame_orf_Score_arr[w]):
                Start = Starts[w]*3 + j +1
                End = temp_stop*3 + j + 3
                NucleotideSeq_Arr = NucleotideSeq_Arr[w]
                CodenLengthArr = CodenLengthArr[w]
    else:
        frame_max_score = 'null'
        Start = 'null'
        End = 'null'
        NucleotideSeq_Arr = 'null'
        CodenLengthArr = 'null'
    return (frame_max_score,Start,End,NucleotideSeq_Arr,CodenLengthArr)
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
            number = float(number) + float(hash_score[temp1])  ########################################## score matrix is very imporant
    mlcds_score.append(number)
    return (mlcds_score)
def mainProcess(inputFile, codonArr, hash_matrix):
    label_Arr, FastA_seq_Arr = Reading(inputFile)
    ######################################################################
    for i in range(len(label_Arr)):
        Label = label_Arr[i]
        Seq = FastA_seq_Arr[i]
        print(i)
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
        ######################################################################six kinds of sequence ######################################################################
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
        ######################################################################
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
        ######################################################################sequence array######################################################################
        Sub_Frame = []
        Seq_len = len(cis_seq_Arr)
        if Seq_len/2 > 61:
            for j in range(0, 6):  ###################################################################### three kinds of open reading frame of each sequence
                TempStr = ''
                if j < 3 :
                    TempStr = InitCodonSeq(j,Seq_len-2,3,cis_seq_Arr)
                    #Total_str_ARR.append(TempStr)
                if 2 < j < 6 :
                    TempStr = InitCodonSeq(j-3,Seq_len-2,3,Re_dou_Arr)
                    #Total_str_ARR.append(TempStr)
                Each_frame_score = CircSimulation(TempStr,hash_matrix)
                Sub_Frame_orf_Score,Sub_Frame_orf,Sub_Start,Sub_End = MaxSubseqPlus(Each_frame_score,TempStr)#########################zxy通过end判断MLCDS是否跨越53端
#                be_len = int(Sub_Start*3 + j)
#                Seq_length = len(Seq_Arr)
#                if be_len > Seq_length:
#                    TempStr = InitCodonSeq(j,Seq_length-2,3,Seq_Arr)
#                    Each_frame_score = CircSimulation(TempStr,hash_matrix)
#                    Sub_Frame_orf_Score,Sub_Frame_orf,Sub_Start,Sub_End = MaxSubseqPlus(Each_frame_score,TempStr)
##################################################################################################################################
                Sub_Frame_orf_seq = Sub_Frame_orf.split(' ')
                Sub_Frame_orf_seq.pop()
                MLCDS_length = len(Sub_Frame_orf_seq) * 3
                CDS_seq = "".join(Sub_Frame_orf_seq)
                MLSCDS_score = CircSimulationmlcds(Sub_Frame_orf,hash_matrix)
                MLSCDS_score_mean = MLSCDS_score[0] / MLCDS_length
                Sub_Frame.append([MLSCDS_score_mean,CDS_seq,MLCDS_length,Sub_Start,Sub_End,j])
                # print(Sub_Frame)####################################################################在这里筛选多个MLCDS
            MLSCDS_score_mean = list(map(lambda x:x[0],Sub_Frame))
            max_index = MLSCDS_score_mean.index(max(MLSCDS_score_mean))
##################################################################################################################
            MLSCDS_score_mean,CDS_seq,MLCDS_length,Sub_Start,Sub_End,j = Sub_Frame[max_index]
            if j <= 2:
                score, number = add_m6a_info(Seq)
                feature_score = get_feature_score(Seq)
            else:
                score, number = add_m6a_info(re_Seq)
                feature_score = get_feature_score(re_Seq)
            lengthbi = MLCDS_length / Seq_len#################################qiqi
            Frame_start_coden_sequence = []
            Frame_stop_coden_sequence = []
            Directory = ''
            if j < 3:
                Directory = '+'
            else:
                Directory = '-'
            CDS_Score.append(MLSCDS_score_mean)####################################################################################################zxy_MLCDS二连密码子得分
            CDS_length.append(MLCDS_length)##########################################################################################zxy_MLCDS长度
            CDS_sequence.append(CDS_seq)######################################################################################################zxy_MLCDS序列
            ML_RL.append(lengthbi)#################################qiqi
            Seq_len = len(Seq_Arr)
            if j < 3 :######################################################################################zxy
                Three_ref_arr = InitCodonSeq(j,Seq_len*3-2,3,Three_ref_arr)######################################################################################zxy
                Four_ref_arr = InitCodonSeq(j,Seq_len*4-2,3,Four_ref_arr)######################################################################################zxy
                Seven_ref_arr = InitCodonSeq(j,Seq_len*7-2,3,Seven_ref_arr)######################################################################################zxy
                Eight_ref_arr = InitCodonSeq(j,Seq_len*8-2,3,Eight_ref_arr)######################################################################################zxy
            if 2 < j < 6 :######################################################################################zxy
                Three_ref_arr = InitCodonSeq(j,Seq_len*3-2,3,Re_Three_ref_arr)######################################################################################zxy
                Four_ref_arr = InitCodonSeq(j,Seq_len*4-2,3,Re_Four_ref_arr)######################################################################################zxy
                Seven_ref_arr = InitCodonSeq(j,Seq_len*7-2,3,Re_Seven_ref_arr)######################################################################################zxy
                Eight_ref_arr = InitCodonSeq(j,Seq_len*8-2,3,Re_Eight_ref_arr)######################################################zxy
            Frame_start_coden_sequence, Frame_stop_coden_sequence = FindStartCodne(Three_ref_arr,Four_ref_arr,Seven_ref_arr,Eight_ref_arr,Sub_Start,Sub_End,j,Seq_len)
#######################################################################candidate sorf information#############################################################
            if len(Frame_start_coden_sequence) >= 1 and len(Frame_stop_coden_sequence) >= 1:
                FrameScore, FrameStart, FrameEnd, nucleotideSeq, CodenSeqLength = Small_orf_score_max(Three_ref_arr,Four_ref_arr,Seven_ref_arr,Eight_ref_arr,Frame_start_coden_sequence,Frame_stop_coden_sequence,hash_matrix,j,Sub_End,Seq_len)
                RnaOrfSequence.append(nucleotideSeq)
                RnaOrfStart.append(FrameStart)
                RnaOrfStop.append(FrameEnd)
                RnaOrfLength.append(CodenSeqLength)
                RnaOrfScore.append(FrameScore)
                RnaOrfFrameIndex.append(j)
                RnaOrfDirectory.append(Directory)
            else:
                RnaOrfSequence.append('null')
                RnaOrfStart.append('null')
                RnaOrfStop.append('null')
                RnaOrfLength.append('null')
                RnaOrfScore.append('null')
                RnaOrfFrameIndex.append(j)
                RnaOrfDirectory.append(Directory)
#######################################################################################
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
                '/home/ZXY/software/PhyloCSF/PhyloCSF 100vertebrates --strategy=fixed --minCodons=30 --frames=3 --removeRefGaps  ' + CDS_path_name + ' |awk \'{print $3}\'').read()
            Conservation_Score_Read = Conservation_Score_Read.replace('\n', '')
            CDS_conservation.append(Conservation_Score_Read)
#######################################################################################################
            RnaOutInformation = ''
            RnaOutInformation = RnaOutInformation + 'Index:' + str(RnaOrfFrameIndex[0]) + '\t'
            RnaOutInformation = RnaOutInformation + 'Directory:' + str(RnaOrfDirectory[0]) + '\t'
            RnaOutInformation = RnaOutInformation + 'Sequence:' + str(RnaOrfSequence[0]) + '\t'
            RnaOutInformation = RnaOutInformation + 'Score:' + str(RnaOrfScore[0]) + '\t'
            RnaOutInformation = RnaOutInformation + 'Start:' + str(RnaOrfStart[0]) + '\t'
            RnaOutInformation = RnaOutInformation + 'End:' + str(RnaOrfStop[0]) + '\t'
            RnaOutInformation = RnaOutInformation + 'CDS_Score:' + str(MLSCDS_score_mean) + '\t'
            RnaOutInformation = RnaOutInformation + 'Length:' + str(RnaOrfLength[0])
            RnaOutInformation = 'CircRNA_ID:' + str(Label) + '\t' + RnaOutInformation + '\n'
            Temp_Out_Information.write(RnaOutInformation)
            ###################################################################### over######################################################################
            CDS_Feature = ''
            CDS_Feature = CDS_Feature + str(CDS_conservation[0]) + ' '
            CDS_Feature = CDS_Feature + str(CDS_Score[0]) + ' '
            CDS_Feature = CDS_Feature + str(CDS_length[0]) + ' '
            CDS_Feature = CDS_Feature + str(ML_RL[0]) + ' '###########################qiqi
            ###################################################################### add m6a information
            CDS_Feature = CDS_Feature + str(number) + ' '
            CDS_Feature = CDS_Feature + str(score) + ' '
            CDS_Feature = CDS_Feature + str(feature_score) + ' '
            CDS_Feature = CDS_Feature + str(CDS_sequence[0])
            CDS_Feature = str(Label) + '|' + str(j) + ' ' + CDS_Feature + ' ' + str(SampleClase) + '\n'
            Temp_Out_Feature.write(CDS_Feature)
            if i == len(label_Arr)-1:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}")
            else:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}\r\n")
        else:
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
#####################################################################################################
            n = 0
            RnaOutInformation = ''
            RnaOutInformation = RnaOutInformation + 'Index:' + str(RnaOrfFrameIndex[0])  + '\t'
            RnaOutInformation = RnaOutInformation + 'Directory:' + str(RnaOrfDirectory[0]) + '\t'
            RnaOutInformation = RnaOutInformation + 'Sequence:' + str(RnaOrfSequence[n]) + '\t'
            RnaOutInformation = RnaOutInformation + 'Score:' + str(RnaOrfScore[n]) + '\t'
            RnaOutInformation = RnaOutInformation + 'Start:' + str(RnaOrfStart[0]) + '\t'
            RnaOutInformation = RnaOutInformation + 'End:' + str(RnaOrfStop[n]) + '\t'
            RnaOutInformation = RnaOutInformation + 'CDS_Score:' + str(MLSCDS_score_mean) + '\t'##################################################################################???是MLCDS_score还是CDS_score
            RnaOutInformation = RnaOutInformation + 'Length:' + str(RnaOrfLength[n])
            RnaOutInformation = 'LncRNA_ID:' + str(Label) + '\t' +  RnaOutInformation + '\n'
            Temp_Out_Information.write(RnaOutInformation)
########################################################################################################
            CDS_Feature = ''
            CDS_Feature = CDS_Feature + str(CDS_conservation[0]) + ' '
            CDS_Feature = CDS_Feature + str(CDS_Score[0]) + ' '
            CDS_Feature = CDS_Feature + str(CDS_length[0]) + ' '
            CDS_Feature = CDS_Feature + str('0') + ' '###########################qiqi
            ########################################## add m6a information
            CDS_Feature = CDS_Feature + str(number) + ' '
            CDS_Feature = CDS_Feature + str(score) + ' '
            CDS_Feature = CDS_Feature + str(feature_score) + ' '
            CDS_Feature = CDS_Feature + str(CDS_sequence[0])
            CDS_Feature = str(Label) + '|' + str('0') + ' ' + CDS_Feature + ' ' + str(SampleClase) + '\n'
            Temp_Out_Feature.write(CDS_Feature)
            if i == len(label_Arr) - 1:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}")
            else:
                Temp_Out_Feature_fasta.write(f"{Label}\r\n{Seq.upper()}\r\n")
######################################################################e:
import uuid
import shutil
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
    parse = optparse.OptionParser()
    # FileName = raw_input('Please enter your a file name: ')
    parse.add_option('-f', '--file', dest='file', action='store', metavar='input files',
                     help='enter your transcript (sequence or gtf)')
    parse.add_option('-o', '--out', dest='outfile', action='store', metavar='output files',
                     help='assign your output file')
    parse.add_option('-s', '--sample', dest='sample', action='store', metavar='sample class', default='1',
                     help='please enter your specified speed ratio')
    parse.add_option('-m', '--model', dest='model', action='store', metavar='model types', default='ve',
                     help='please enter your specified classification model')
    # parse.add_option('-g','--gtf',dest='gtf',action='store_true',metavar='gtf file name',help='please enter your gtf files')
    # parse.add_option('-d','--directory',dest='directory',action='store',metavar='',help='if your input file is gtf type please enter RefGenome directory')
    (options, args) = parse.parse_args()
    inPutFileName = options.file
    outPutFileName = options.outfile
    SampleClase = int(options.sample)
    ClassModel = options.model
    # over######################################################################
    curPath = cur_file_dir()
    MatrixPath = "/home/ZXY/peptide/circRNA/score_matrix/cir_score_matrix"
#    MatrixPath += "/score_matrix/cir_score_matrix"
    inMatrix = open(MatrixPath)
    Matrix = inMatrix.read()
    inMatrix.close()
    Circ_Dir = outPutFileName + '_Tmp_Dir'
    subprocess.call('mkdir ' + Circ_Dir + '', shell=True)
    Out_Information = Circ_Dir + '/Information'
    Temp_Out_Information = open(Out_Information, 'w')
    Out_Feature = Circ_Dir + '/Feature'
    Temp_Out_Feature = open(Out_Feature, 'w')
    Out_Feature_fasta = Out_Feature + "_fasta.fa"
    Temp_Out_Feature_fasta = open(Out_Feature_fasta, 'w')
    CDS_path = Circ_Dir + '/C_score_write'
    # Conservation_Score_Read = Circ_Dir + '/C_score_Read'
    # Temp_Con_Score_Read = open(Conservation_Score_Read,'w')
    # over######################################################################
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
    # over######################################################################
    inFiles = open(inPutFileName)
    inFilesArr = inFiles.read()
    inFileNum = inFilesArr.split('\n')
    inFileLen = len(inFileNum) - 1
    inFiles.close()
    sequence_Arr = inFilesArr.split('\n')  ###### input data
    # sLen = len(sequence_Arr) - 1
    # del sequence_Arr[sLen]
    ARRAY = TwoLineFasta(sequence_Arr)  ####### transform multiple line to Two line sequence
    Label_Array, FastA_Seq_Array = Tran_checkSeq(ARRAY)
    inFileLength = len(Label_Array)
    TOT_STRING = []
    # over######################################################################label modification######################################################################
    for i in range(len(Label_Array)):
        tmp_label_one = Label_Array[i]
        tmp_label = tmp_label_one.replace('\r', '')
        tmp_seq = FastA_Seq_Array[i]
        Temp_Seq = tmp_seq.replace('\r', '')
        TOT_STRING.append(tmp_label)
        TOT_STRING.append(Temp_Seq)
    ######################################################################input sequence array######################################################################
    mainProcess(TOT_STRING, Alphabet, Matrix_hash)
    # over######################################################################
    Temp_Out_Information.close()
    Temp_Out_Feature.close()
    Temp_Out_Feature_fasta.close()
    # Temp_Con_Score_Write.close()
    Train_Info = open(Out_Information)
    Train_Info_Seq = Train_Info.read()
    Train_Info_Arr = Train_Info_Seq.split('\n')
    Train_Feature = open(Out_Feature)
    Train_Feature_Seq = Train_Feature.read()
    Train_Feature_Arr = Train_Feature_Seq.split('\n')
    Pre_mian('%s/Feature' % Circ_Dir, Out_Information,ClassModel)