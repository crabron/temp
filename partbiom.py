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
import time
import collections


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

def benchmark(func):
    def wrapper(*args, **kwargs):
        t = time.clock()
        res = func(*args, **kwargs)
        print(func.__name__, time.clock() - t)
        return res
    return wrapper

@benchmark
def del_singltones(inp):
    otu_table = biom.load_table(inp)
    a = [i for i in otu_table.iter_data(axis="observation")]
    b = [i for i in a if sum(i)>5]

@benchmark
def ratiometr(ftable, ilist):  
    left = [i for i in ftable.columns.values if i.startswith(ilist[0])]
    right = [i for i in ftable.columns.values if i.startswith(ilist[1])]
    sum1 = ftable[left].sum().sum().sum()
    sum2 = ftable[right].sum().sum().sum()
    rtableleft = ftable[left].apply(lambda row: row/sum1)
    rtableright = ftable[right].apply(lambda row: row/sum2)
    table = pandas.concat([rtableleft, rtableright, ftable.taxonomy], axis=1, join='inner')
    return left, right, table, sum1, sum2

@benchmark
def ratfilt(left, right, table):
    table['mean1'] = table[left].apply(lambda row: sum(row)/len(row),axis=1)
    table['mean2'] = table[right].apply(lambda row: sum(row)/len(row),axis=1)
    ftable = table.ix[(table['mean1'] >= 0.005) | (table['mean2'] >= 0.005)]
    return ftable


@benchmark
def man_y(rttable,left, right):
    def many(x,y):
        f = man(x,y, use_continuity=True, alternative=None)[1]
        return f
    
    rttable['pvalue'] = rttable.apply(lambda x: many(x[left],x[right]),axis=1)
    return rttable

@benchmark
def filter_otu(inp, idfs):
    '''
    Filter biom table by sample Id from mapping txt file.
    '''
    otu_table = biom.load_table(inp)
    slist = otu_table.ids(axis='sample')
    ll = [i for i in slist if i.startswith(idfs[0]) or i.startswith(idfs[1])]
    new_table = otu_table.filter(ll,axis='sample', inplace=False)
    return new_table

@benchmark
def log(out, table,ilist):
    idfs = '_'.join(ilist)
    summ = table.count()[0]
    t_unic1 = table[(table.mean1 > 0) & (table.mean2 == 0)].count()[0]
    t_unic2 = table[(table.mean2 > 0) & (table.mean1 == 0)].count()[0]
    table_gen = table.ix[(table.mean1 > 0) & (table.mean2 > 0)]
    general = table_gen.count()[0]
    increase_in_1 = table_gen[(table.mean1 > table.mean2) & table.pvalue<= 0.05].count()[0]
    decrease_in_1 = table_gen[(table.mean1 < table.mean2) & table.pvalue<= 0.05].count()[0]
    same = table_gen[table.pvalue > 0.05].count()[0]
    a = [summ, t_unic1, t_unic2, general ,increase_in_1, decrease_in_1, same]
    out.loc[idfs] = a
    return out, table_gen

@benchmark
def out(ilist, table, sum1, sum2, left, right, table_gen):
    a = "_".join(ilist)
    writer = pandas.ExcelWriter('script/{}.xlsx'.format(a))
    t_unic1 = table.ix[(table.mean1 > 0) & (table.mean2 == 0)]
    df2 = t_unic1[['taxonomy','mean1','mean2']]
    t_unic2 = table.ix[(table.mean2 > 0) & (table.mean1 == 0)]
    df3 = t_unic2[['taxonomy','mean1','mean2']]
    table_gen_change = table_gen.ix[table.pvalue<= 0.05]
    table_gen_change['proportion'] = table_gen_change.apply(lambda x: x['mean1']/x['mean2'],axis=1)
    df1 = table_gen_change[['taxonomy','mean1','mean2','proportion']]
    df1.to_excel(writer,'general_with_proportion')
    df2.to_excel(writer,'unic_for_1')
    df3.to_excel(writer,'unic_for_2')
    writer.save()
    print()
    # table_change = table.ix[table.pvalue <= 0.05]
    # table_change['pseudo_mean1'] = table_change[left].apply(lambda row: row/sum1)
    # d_out = pandas.DataFrame(columns=('proportion' ,'prop_in_1',  'prop_in_2')) 
    # ttable_change['proportion'] = table_change.apply(lambda x: x['mean1']/['mean2'],axis=1)
    return

@benchmark
def transform(table):
    out = table.to_dataframe()
    return out

@benchmark
def alltransform(table):
    ftable = biom.load_table(table)
    d = []
    for i,e,w in ftable.iter(dense=True, axis='observation'):
        for a in w.values():
            # ar = a[::-1]
            vall = "|".join(a)
            d.append(vall)
    out = ftable.to_dataframe()
    out['taxonomy'] = d
    return out

@benchmark
def pfilter(table, ilist):
    ids = [i for i in table.columns.values if i.startswith(ilist[0]) or i.startswith(ilist[1])]
    ids.append('taxonomy')
    ftable = table[ids]
    return(ftable)

def main(inp, ID):
    try:
        if not os.path.exists("script"):
            os.makedirs("script")
        d_out = pandas.DataFrame(columns=('summ' ,'unic_for_1',  'unic_for_2', 'general', 'increase_in_1', 'decrease_in_1', 'same'))
        with open(ID, "r") as idfile:
            idf = idfile.readlines()
            lenidf = len(idf)
            if lenidf >= 2:
                table = alltransform(inp)
            for i in idf:
                nonn_line = i.rstrip()
                ilist = nonn_line.split("\t")
                if lenidf <= 1:
                    ftable = filter_otu(inp, ilist)
                    datatable = transform(ftable)
                else:
                    datatable = pfilter(table, ilist)
                left, right, rttable, sum1, sum2 = ratiometr(datatable, ilist)
                frttable = ratfilt(left, right, rttable)
                frttable = man_y(frttable,left, right)
                global table_out
                table_out, table_gen = log(d_out, frttable, ilist)
                out(ilist, frttable, sum1, sum2, left, right, table_gen)
        with open("script/log.csv", "w") as l:
            table_out.to_csv( l , sep='\t')
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
