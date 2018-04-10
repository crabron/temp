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
            w = [ftable.get_value_by_ids(q, a) for a in i]
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

def fsometests(sumotudict, ftable, iddict, sumalldict, otu_list, whattest="chi2"):

    '''
    This function consists of statistical tests, that calculate p-value for our data.
    This is my  shiny shit castle of crap, which really lacks some order.
    Return orderdicts with pairs like  otu:p-value: 
    fisher/chi2(if sum from some sample more than 5)[0]
    ttest [1]
    kruskal[2]
    '''

    otuwsa_od = OrderedDict()
    for i in otu_list:
        dpartbef = [a for a in iddict.values()]
        dpart1 = [ftable.get_value_by_ids(i, a) for a in dpartbef[0]]
        dpart2 = [ftable.get_value_by_ids(i, a) for a in dpartbef[1]]
        otuwsa_od.update({i:[dpart1,dpart2]})



    otherdict = OrderedDict()
    sumalldict = sumalldict.values()
    for i in sumotudict.items():
        other = map(operator.sub, sumalldict , i[1])
        otherdict.update({i[0]:other})
    leno = len(otu_list)
   
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
            if a == 0: 
                a == 0.001

            if b == 0:
                b == 0.001
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

    pdict = fsometests(f_sumotudict, ftable_ab, f_iddict, sumalldict, f_otu_list, whattest)
    ftable_ab_pv = ftable_del_pval(pdict, ftable_ab)

    f_iddict = fiddict(ftable_ab_pv)
    f_otu_list = ftable_ab_pv.ids(axis='observation')
    f_sumotudict = fsumotudict(ftable_ab_pv, f_otu_list, f_iddict)


    fdiv, part = fbiom_differ(f_sumotudict, sumalldict)
    ftable_ab_pa = fdelete_minors_part(ftable_ab_pv, part, part_in)

    a_iddict = fiddict(ftable_ab_pa)
    a_otu_list = ftable_ab_pa.ids(axis='observation')
    a_sumotudict = fsumotudict(ftable_ab_pa, a_otu_list, a_iddict)
    # a_sumotudict_pseudo = pseudoc(a_sumotudict)
    div, part = fbiom_differ(a_sumotudict, sumalldict)
    return div, part
    
# if __name__ == "__main__":
#     main(inp, ID, part_in)
    


def USELESS():
    '''mainflow'''

    otu_table = biom.load_table(inp) 
    ftable = filter_otu(otu_table, idfs) #filtered biom table
    sample_list = ftable.ids(axis='sample') #all sample ids
    otu_list = ftable.ids(axis='observation') #all obs ids(otus)
    f_iddict = iddict(sample_list)  #groupped otu sample: key - master, value - master replication
    otuwsa_od = otsa_od(otu_list,f_iddict) 
    sumotudict = summotu_od(otu_list,f_iddict)
    sumalldict = summall_od(otu_list,f_iddict, sumotudict)


    ''' part with hop, hop pseudozeros and minor otu eliminating'''

    sys.stdout.write('\r')
    sys.stdout.write("Doing divisions and subtrations for greater good.")
    sys.stdout.flush()


    
    fddiv, fdsub, part = fbiom_differ(sumotudict, sumalldict)
    table_del_ab = table_del_ab(ftable)
    pdict, tdict, kdict, mdict = sometests(sumotudict,otuwsa_od)
    delete_minors_part(table, part, part_in)
    fd = table_del_pval(p_dict, table_del) 
    s_otu_list = ftable_pval.ids(axis='observation')
    fsample_list = ftable_pval.ids(axis='sample')
    s_otu_list = ftable_pval.ids(axis='observation')
    s_iddict = iddict(fsample_list) 
    otuwsa_od = otsa_od(s_otu_list,s_iddict) 
    f_sumotudict = summotu_od(s_otu_list,s_iddict) 
    fbd = fbiom_differ(sumotudict, f_sumotudict, sumalldict)

    sumalldict = summall_od(otu_list,f_iddict, sumotudict)
    pvalues = sometests(sumotudict,otuwsa_od) 
    p_dict = pvalues[0] 
    fd = table_del_pval(p_dict, ftable) 


    '''start doing monkey job again for new filtered OTU table'''

    sys.stdout.write('\r')
    sys.stdout.write("Doing divisions and subtrations for greater good.")
    sys.stdout.flush()

    fsample_list = fd.ids(axis='sample')
    s_otu_list = fd.ids(axis='observation')
    s_iddict = iddict(fsample_list) 
    otuwsa_od = otsa_od(s_otu_list,s_iddict) 
    f_sumotudict = summotu_od(s_otu_list,s_iddict) 
    fbd = fbiom_differ(sumotudict, f_sumotudict, sumalldict)

    print("otu", "fisher or chi", "t-test","krus", sep="\t", end="\n", file=results)
    for i in otu_list:
        pd = pvalues[0]
        td = pvalues[1]
        kd = pvalues[2]
        fp = pd.get(i)
        tp = td.get(i)
        kp = kd.get(i)
        otusam1 = otuwsa_od.get(i)[0]
        otusam2 = otuwsa_od.get(i)[1]
        print(i, fp, tp, kp, sep="\t", end="\n", file=pval)
        '''out'''




# if not os.path.exists("script"):
#     os.makedirs("script")

if not os.path.exists("script"):
    os.makedirs("script")
with open(ID, "r") as idfile:
    idf = idfile.readlines()
    for i in idf:
        div, part = main(inp, i, part_in)
        d = i.rstrip()
        b = d.split("\t")
        a = "_".join(b)
        with open("script/subdel_{}.txt".format(a), "a") as out:
            print("otu", "division","part_1","part_2",whattest, sep="\t", end="\n", file=out)
            for i in part.keys():
                fdiv = div.get(i)
                p1 = part.get(i)[0]
                p2 = part.get(i)[1]
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

