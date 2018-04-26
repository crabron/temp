#! /usr/bin/env python

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu as man
import biom
import pandas
import h5py
import re
import itertools
import sys
import os
import shutil
import traceback

inp = argv[1]
if argv[2] == 'all':
    otu_table = biom.load_table(inp)
    sample_list = otu_table.ids(axis='sample')
    atata = []
    for a in sample_list:
        i = a.split('.')
        ca = i[:-1]
        caa = ".".join(ca)
        atata.append(caa)
    at_set = set(atata)
    w=[i for i in itertools.combinations(at_set, 2)]
    with open('temp2.txt', 'a') as ID:
        for i in w:
            ID.write(i[0] + '\t' + i[1]+ '\n')
    ID = 'temp2.txt'
else:
    ID = argv[2]

def del_singltones(inp):
    otu_table = biom.load_table(inp)
    a = [i for i in otu_table.iter_data(axis="observation")]
    b = [i for i in a if sum(i)>5]

def ratiometr(ftable, ilist):  
    left = [i for i in ftable.columns.values if i.startswith(ilist[0])]
    right = [i for i in ftable.columns.values if i.startswith(ilist[1])]
    tableleft = ftable.groupby(left)
    tableright = ftable.groupby(right)
    sum1 = tableleft.sum().sum().sum()
    sum2 = tableright.sum().sum().sum()
    rtableleft = ftable[left].apply(lambda row: row/sum1)
    rtableright = ftable[right].apply(lambda row: row/sum2)
    table = pandas.concat([rtableleft, rtableright], axis=1, join='inner')
    return left, right, table

def ratfilt(left, right, table):
    table['mean1'] = table[left].apply(lambda row: sum(row)/len(row),axis=1)
    table['mean2'] = table[right].apply(lambda row: sum(row)/len(row),axis=1)
    table['tresh_ratio'] = table['mean1' or 'mean2'] >=0.005
    return table


def man_y(rttable,left, right):
    def many(x,y,t):
        if t == True:
            f = man(x,y, use_continuity=True, alternative=None)[1]
            return f
    rttable['pvalue'] = rttable.apply(lambda x: many(x[left],x[right],x['tresh_ratio']),axis=1)
    whochange = rttable['pvalue'] <= 0.05
    whosame = rttable['pvalue'] > 0.05
    print(rttable)
    return rttable, whochange, whosame

def filter_otu(inp, idfs):
    '''
    Filter biom table by sample Id from mapping txt file.
    '''
    otu_table = biom.load_table(inp)
    al_l =[line.rstrip() for line in idfs.split("\t")]
    al_ll = []
    for a in al_l:
        if re.search(r'\D\Z',a):
            sample_list = otu_table.ids(axis='sample')
            for i in sample_list: 
                if i.startswith(a):
                    al_ll.append(i)
        else:
            al_ll.append(a)
    
    new_table = otu_table.filter(al_ll,axis='sample', inplace=False)
    return new_table

def out(table):
    with open("script/log.txt", "a") as log:
        




        print(idfs, otusum , unic1 ,unic2,  lenlconst, lenlmore_1, lenlless_1, lenlmore_2, lenlless_2, sep="\t", end="\n", file=log)
    return


def main(inp, ID):
    try:
        if not os.path.exists("script"):
            os.makedirs("script")
        with open("script/log.txt", "a") as l:
            print("name_for_pair", "summ", "unic_for_1", "unic_for_2", "general", "increase_in_1","decrease_in_1","increase_in_2","decrease_in_2", sep="\t", end="\n", file=l)
        with open(ID, "r") as idfile:
            idf = idfile.readlines()
            for i in idf:
                nonn_line = i.rstrip()
                ilist = nonn_line.split("\t")
                ftable = filter_otu(inp, i)
                datatable = ftable.to_dataframe()
                left, right, rttable = ratiometr(datatable, ilist)
                ratfiltuu(left, right, rttable)
                rttable, whochange, whosame = man_y(rttable,left, right)
                out(rttable, whochange, whosame)
        if os.path.isfile("temp2.txt"):   
            os.remove('temp2.txt')

    except Exception:
        if os.path.isfile("temp2.txt"):   
            os.remove('temp2.txt')
        if os.path.exists("script"):
            shutil.rmtree('script', ignore_errors=True)
        print(traceback.format_exc())

if __name__ == "__main__":
    main(inp, ID)
