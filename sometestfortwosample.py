#! /usr/bin/env python

"""
Work only for couples. Ignor groups. Relise of preveous code, but more kinda introvertic.
"""

from __future__ import division, print_function
from scipy.stats import fisher_exact as fisher
from scipy.stats import chi2_contingency as chisq
from scipy.stats import ttest_ind as ttest
from scipy.stats import mannwhitneyu as man
from scipy.stats import ranksums as wilc
from scipy.stats import kruskal as krus
from scipy import stats
import getopt
import biom
import re
from collections import OrderedDict
import timeit
import numpy as np
import itertools as it
import operator
import argparse
from time import sleep
import sys
from sys import argv
import pandas
import math
import os
import matplotlib.pyplot as plt

inp = argv[1]
ID = argv[2]
part_in = argv[3]
whattest = argv[4]

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

def fiddict(ftable):
    '''
    Ordered dict with groupped otu sample: key - master, value - master replication
    '''
    samplesame_list = []
    sample_list = ftable.ids(axis='sample')
    for i in sample_list: 
        q = re.match(r'\D+', i)
        samplesame_list.append(q.group())
    sampledd_list = list(OrderedDict.fromkeys(samplesame_list))
    sorted_sample = []
    for i in sampledd_list:
        la = []
        for ia in sample_list:
            if ia.startswith(i):
                la.append(ia)
        sorted_sample.append(la)
    key_id = [re.split(r'.\d',i[0])[0] for i in sorted_sample]
    zipl = zip(key_id,sorted_sample)
    iddict = OrderedDict(zipl)
    return iddict

def fsumotudict(ftable, otu_list, iddict):
    '''
    Create OrderedDict odject
    key = OTU Id
    value = list of sum for all repeats in samle order
    '''
    sumotudict = OrderedDict()
    for q in otu_list:
        suml = []
        for i in iddict.values():
            w = []
            for a in i:
                wa = ftable.get_value_by_ids(q, a)
                w.append(wa)
                # print(q, wa, w)

            # w = [(ftable.get_value_by_ids(q, a)) for a in i]
            sum_otu = sum(w)
            suml.append(sum_otu)
        sumotudict.update({q:suml})
    return sumotudict

def fsumalldict(otu_list,iddict, sumotudict):
    '''
    Function create OrderedDict object.
    key = sample Id
    value = sum of all otu for this sample
    '''
    sumalldict = OrderedDict()
    leng = len(otu_list)
    lenglist = range(0, leng)
    keyids = iddict.keys()
    zipl = zip(keyids, lenglist)
    for w in zipl:
        a = w[1]
        r = sumotudict.values()
        y =[t[a] for t in r]
        sum_y = sum(y)
        sumalldict.update({w[0]:sum_y})
    return sumalldict

def fbiom_differ(f_sumotudict, sumalldict):
    '''counts part, division parts from each other'''
    div = OrderedDict()
    part = OrderedDict()
    for d in f_sumotudict.items():
        key = d[0]
        summ_i = sumalldict.values()
        summ_1 = summ_i[0]
        summ_2 = summ_i[1]
        val_w = d[1]
        a = val_w[0]
        b = val_w[1]
        if a == 0:
            part_1 = 1/summ_1
        else:
            part_1 = a/summ_1

        if b == 0:
            part_2 = 1/summ_2
        else:
            part_2 = b/summ_2
        part.update({key:[part_1, part_2]})
    for d in part.items():
        key = d[0]
        w = d[1]
        if w[0] == 0:
            fdiv = 1 / w[1]
        elif w[1] == 0:
            fdiv = w[0] / 1
        else:
            fdiv = w[0] / w[1]
        div.update({key:fdiv})
    # fdsub = OrderedDict()
    # for d in part.items():
    #     key = d[0]
    #     w = d[1]
    #     sub = w[0] - w[1]
    #     fdsub.update({key:sub})
    return div, part

def fdelete_minors_ab(sumotudict, ftable): 
    otu_del = [ d[0] for d in sumotudict.items() if sumotudict.get(d[0])[0] >= 6 or sumotudict.get(d[0])[1]>=6 ]
    ftable_ab = ftable.filter(otu_del,axis='observation', inplace=False)
    return ftable_ab

def fdelete_minors_part(ftable, part, part_in):
    a = float(part_in)
    otu_del = [ d[0] for d in part.items() if part.get(d[0])[0] >= a or part.get(d[0])[1]>= a]
    ftable_ab = ftable.filter(otu_del, axis='observation', inplace=False) 
    return ftable_ab

def fsometests(sumotudict, ftable, iddict, sumalldict, otu_list, whattest="man-y"):

    '''
    This function consists of statistical tests, that calculate p-value for our data.
    This is my  shiny shit castle of crap, which really lacks some order.
    Return orderdicts with pairs like  otu:p-value: 
    fisher/chi2(if sum from some sample more than 5)[0]
    ttest [1]
    kruskal[2]
    '''
    
    otherdict = OrderedDict()
    sumalldict = sumalldict.values()
    for i in sumotudict.items():
        other = map(operator.sub, sumalldict , i[1])
        otherdict.update({i[0]:other})
    leno = len(otu_list)

    
    otuwsa_od = OrderedDict()
    for i in otu_list:
        dpartbef = [a for a in iddict.values()]
        summm_1 = sumalldict[0]
        summm_2 = sumalldict[1]
        dpart1_l = [ftable.get_value_by_ids(i, a) for a in dpartbef[0]]
        dpart2_l = [ftable.get_value_by_ids(i, a) for a in dpartbef[1]]
        dpart1 = [a/summm_1 for a in dpart1_l]
        dpart2 = [a/summm_2 for a in dpart2_l]
        otuwsa_od.update({i:[dpart1,dpart2]})
   
    if whattest == "chi2":
        fj=0
        pdict = OrderedDict()
        for i in sumotudict.items():
            fj+=1
            sum = i[1]
            other = otherdict.get(i[0])
            table = np.array([sum,
                            other])
            if sum[1] <= 5 or sum[0] <= 5: 
                p = fisher(table)[1]
                pdict.update({i[0]:p})
            else:
                p = chisq(table, lambda_="log-likelihood")[1]
                pdict.update({i[0]:p})
            
            sys.stdout.write('\r')
            sys.stdout.write("fisher {}/{}".format(fj,leno))
            sys.stdout.flush()
    
    elif whattest=="ttest":

        pdict = OrderedDict()
        j=0
        for i in otu_list:
            j+=1
            a = otuwsa_od.get(i)[0]
            b = otuwsa_od.get(i)[1]
            p = ttest(a,b)[1]
            if p != p:
                p = 1
            pdict.update({i:p})

            sys.stdout.write('\r')
            sys.stdout.write("ttest {}/{}".format(j,leno))
            sys.stdout.flush()
        
    elif whattest=="man-y":

        pdict = OrderedDict()
        mj=0
        for i in otu_list:
            mj+=1
            a = otuwsa_od.get(i)[0]
            b = otuwsa_od.get(i)[1]
            p = man(a,b)[1]
            pdict.update({i:p})
            sys.stdout.write('\r')
            sys.stdout.write("man-y {}/{}".format(mj,leno))
            sys.stdout.flush()
    
    elif whattest=="kruscal":

        pdict = OrderedDict()
        kw=0
        for i in otu_list:
            kw+=1
            a = otuwsa_od.get(i)[0]
            b = otuwsa_od.get(i)[1]
            p = wilc(a,b)[1]
            pdict.update({i:p})
            sys.stdout.write('\r')
            sys.stdout.write("kruskal {}/{}".format(kw,leno))
            sys.stdout.flush()
        
        sys.stdout.write('\r')
        
    return pdict

def ftable_del_pval(pdict, new_table):
    fdd = [d[0] for d in pdict.items() if d[1] <= 0.05]
    ftable_pval = new_table.filter(fdd, axis='observation', inplace=False)
    return ftable_pval

def pseudoc(sumotudict):
    ps = OrderedDict()
    for i in sumotudict.items():
        key = i[0]
        val = i[1]
        val_1 = []
        val_2 = []
        if val[0] == 0:
            val_1.append(1)
        else:
            val_1.append(val[0])
        if val[1] == 0:
            val_2.append(1)
        else:
            val_2.append(val[1])
        v = zip(val_1,val_2)
        ps.update({key:v})
    return ps

def log(pdict, div, ftable_ab_pa, idfs, a_sumotudict):
    with open("script/log.txt", "a") as log:

        otus = ftable_ab_pa.ids(axis='sample')
        vs = len(otus)

        lunic1 = [a for a in a_sumotudict.items() if a[1][0] != 0 and a[1][1] == 0] 
        lunic2 = [a for a in a_sumotudict.items() if a[1][1] != 0 and a[1][0] == 0] 
        lenunic1 = dict(lunic1)
        lenunic2 = dict(lunic2)
        unic1 = len(lenunic1)
        unic2 = len(lenunic2)

        lenunic = lenunic1.copy()   
        lenunic.update(lenunic2)  
        nu_sumo = [a for a in a_sumotudict.items() if a[0] not in lenunic.keys()]
        nu_sumotudict = dict(nu_sumo)

        llconst = [a[0] for a in pdict.items() if a[1] >= 0.05]
        lconst = [ a for a in llconst if a in nu_sumotudict.keys()]
        lenlconst = len(lconst)

        llmore = [a[0] for a in div.items() if a[1] >=1]
        lmore = [ a for a in llmore if a in nu_sumotudict.keys()]
        lenlmore = len(lmore)

        lless = [a[0] for a in div.items() if a[1] <=1]
        less = [ a for a in lless if a in nu_sumotudict.keys()]
        lenlless = len(less)

        print(idfs, unic1 ,unic2,  lenlconst, lenlmore, lenlless, sep="\t", end="\n", file=log)
    return

def main(inp, ID, part_in):

    idfs = ID

    sys.stdout.write('\r')
    sys.stdout.write("Prepare biom table for analysis")
    sys.stdout.flush()


    ftable = filter_otu(inp, idfs) #filtered biom table
    iddict = fiddict(ftable)
    otu_list = ftable.ids(axis='observation')
    sumotudict = fsumotudict(ftable, otu_list, iddict)
    sumalldict = fsumalldict(otu_list, iddict, sumotudict)
    ftable_ab = fdelete_minors_ab(sumotudict, ftable)

    f_iddict = fiddict(ftable_ab)
    f_otu_list = ftable_ab.ids(axis='observation')
    f_sumotudict = fsumotudict(ftable_ab, f_otu_list, f_iddict)

    fdiv, part = fbiom_differ(f_sumotudict, sumalldict)
    ftable_ab_pa = fdelete_minors_part(ftable_ab, part, part_in)

    fa_iddict = fiddict(ftable_ab_pa)
    fa_otu_list = ftable_ab_pa.ids(axis='observation')
    fa_sumotudict = fsumotudict(ftable_ab_pa, fa_otu_list, fa_iddict)

    pdict = fsometests(fa_sumotudict, ftable_ab_pa, fa_iddict, sumalldict, fa_otu_list, whattest)
    ftable_ab_pv = ftable_del_pval(pdict, ftable_ab_pa)


    a_iddict = fiddict(ftable_ab_pv)
    a_otu_list = ftable_ab_pv.ids(axis='observation')
    a_sumotudict = fsumotudict(ftable_ab_pv, a_otu_list, a_iddict)
    # a_sumotudict_pseudo = pseudoc(a_sumotudict)
    div, part = fbiom_differ(a_sumotudict, sumalldict)
    return div, part, pdict, ftable_ab_pa, fa_sumotudict
    
# if __name__ == "__main__":
#     main(inp, ID, part_in)
    


def fancies(div, part):
    a = div.values()
    b = [ i[0] for i in part.values()]
    arr = np.array([a]+[b])
    df = arr.DataFrame.plot.scatter(columns=['a', 'b'])
    df.plot.scatter(x='a', y='b');
    print(arr)

# if not os.path.exists("script"):
#     os.makedirs("script")
if not os.path.exists("script"):
    os.makedirs("script")
with open("script/log.txt", "a") as l:
    print("name", "unic_1", "unic_2", "allsame", "more_in1", "less_in1", sep="\t", end="\n", file=l)
with open(ID, "r") as idfile:
    idf = idfile.readlines()
    for i in idf:
        div, part, pdict, ftable_ab_pa, fa_sumotudict = main(inp, i, part_in)
        d = i.rstrip()
        b = d.split("\t")
        a = "_".join(b)
        log(pdict, div, ftable_ab_pa, a, fa_sumotudict)
        with open("script/subdel_{}.txt".format(a), "w") as out:
            print("otu", "division","part_1","part_2",whattest, sep="\t", end="\n", file=out)
            for i in part.keys():
                ffdiv = div.get(i)
                p1 = part.get(i)[0]
                p2 = part.get(i)[1]
                if ffdiv <= 1:
                    fdiv = 1/ffdiv
                else:
                    fdiv = ffdiv
                    
                print(i, fdiv, p1, p2, sep="\t", end="\n", file=out)

    # print("fisher, chi ={}".format(qf),"ttest ={}".format(qt),  sep="\t", end="\n", file=log )

sys.stdout.write('\r')
sys.stdout.write("\r")
sys.stdout.flush()

    # hmn = [otu_table.get_value_by_ids('denovo228',i) for i in u]
    # i = otu_table.iter_data()
        # ind = sumotudict.keys().index(i)
        # indm = ind + 1 in enumerate(zip(otu_list[:-1], otu_list[1:]))
        # sum2 = sumotudict.values()[indm]
        

''' Some fancy plots. Try to realize this part on pandas libraries '''

