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




def filter_otu(otu_table, idfs):
    '''
    Filter biom table by sample Id from mapping txt file.
    '''
    al = idfs.read()
    al_l =[line.rstrip() for line in al.split("\t")]
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


def group_sample(sample_list):
    '''
    Return list of sample ids, with groupped in one list repetitions.
    '''
    samplesame_list = []
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
    return sorted_sample


def summotu_od(otu_list,iddict):
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

def summall_od(otu_list,iddict, sumotudict):
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

def otsa_od(otu_list,iddict):
    '''
    Ordered dict with a pair lists in values for tested samples.
    '''
    otsa_od = OrderedDict()
    for i in otu_list:
        dpartbef = [a for a in iddict.values()]
        dpart1 = [ftable.get_value_by_ids(i, a) for a in dpartbef[0]]
        dpart2 = [ftable.get_value_by_ids(i, a) for a in dpartbef[1]]
        otsa_od.update({i:[dpart1,dpart2]})
    return otsa_od

def sometests(sumotudict,otuwsa_od):

    '''
    This function consists of statistical tests, that calculate p-value for our data.
    This is my  shiny shit castle of crap, which really lacks some order.
    Return orderdicts with pairs like  otu:p-value: 
    fisher/chi2(if sum from some sample more than 5)[0]
    ttest [1]
    kruskal[2]
    '''

    otherdict = OrderedDict()
    summall = sumalldict.values()
    for i in sumotudict.items():
        other = map(operator.sub, summall , i[1])
        otherdict.update({i[0]:other})

    # # bonff = 
    leno = len(otu_list)
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

    tdict = OrderedDict()
    j=0
    for i in otu_list:
        j+=1
        a = otuwsa_od.get(i)[0]
        b = otuwsa_od.get(i)[1]
        p = ttest(a,b)[1]
        if p != p:
            p = 1
        tdict.update({i:p})

        sys.stdout.write('\r')
        sys.stdout.write("ttest {}/{}".format(j,leno))
        sys.stdout.flush()

    # mdict = OrderedDict()
    # mj=0
    # for i in otu_list:
    #     mj+=1
    #     a = otuwsa_od.get(i)[0]
    #     b = otuwsa_od.get(i)[1]
    #     if np.sum(a) == 0 or np.sum(b) == 0:
    #         p = "nan"
    #     else:
    #         p = man(a,b)[1]
    #     mdict.update({i:p})
    #     sys.stdout.write('\r')
    #     sys.stdout.write("man-y {}/{}".format(mj,leno))
    #     sys.stdout.flush()

    kdict = OrderedDict()
    kw=0
    for i in otu_list:
        kw+=1
        a = otuwsa_od.get(i)[0]
        b = otuwsa_od.get(i)[1]
        p = wilc(a,b)[1]
        kdict.update({i:p})
        sys.stdout.write('\r')
        sys.stdout.write("kruskal {}/{}".format(kw,leno))
        sys.stdout.flush()
    
    sys.stdout.write('\r')
    
    return pdict, tdict, kdict

def biom_rel_filter(pdict, new_table):
    fdd = [d[0] for d in pdict.items() if d[1] <= 0.05]
    fd = new_table.filter(fdd, axis='observation', inplace=False)
    return fd

# def fbiom_differ(fd,iddict,sumotudict):
#     sumo = np.sum(sumotudict.values())
#     otu = sumotudict.keys()
#     print(iddict)
#     print(otu)
#     fbd = [values[1] for values, id, metadata in fd.iter()]
#     # part_fd = [ d/sumo for d in fd.get_value_by_ids()]
#     return fbd
    


inp = "nice_little_table.biom"
idfs = open("IDtestlist.txt", "r")
otu_table = biom.load_table(inp) 
ftable = filter_otu(otu_table, idfs) #filtered biom table
sample_list = ftable.ids(axis='sample') #all sample ids
otu_list = ftable.ids(axis='observation') #all obs ids(otus)
gs = group_sample(sample_list)  #groupped
key_id = [re.split(r'.\d',i[0])[0] for i in gs]
zipl = zip(key_id,gs)
iddict = OrderedDict(zipl)
otuwsa_od = otsa_od(otu_list,iddict)
sumotudict = summotu_od(otu_list,iddict)
sumalldict = summall_od(otu_list,iddict, sumotudict)
pvalues = sometests(sumotudict,otuwsa_od)
p_dict = pvalues[0]
fd = biom_rel_filter(p_dict, ftable)
# fb_d = fbiom_differ(fd, iddict, sumotudict)
print(gs)
print(zipl)

# start doing monkey job again



# with open("pval.txt", "w") as pval:
#     print("otu", "fisher or chi", "t-test","krus", sep="\t", end="\n", file=results)
#     for i in otu_list:
#         pd = pvalues[0]
#         td = pvalues[1]
#         kd = pvalues[2]
#         fp = pd.get(i)
#         tp = td.get(i)
#         kp = kd.get(i)
#         otusam1 = otuwsa_od.get(i)[0]
#         otusam2 = otuwsa_od.get(i)[1]
#         print(i, fp, tp, kp, sep="\t", end="\n", file=pval)













# with open("log.txt", "w") as log:
#     qf=0
#     qt=0
#     for i in otu_list:
#         pdict = pvalues[0]

#         fp = pdict.get(i)
#         tp = tdict.get(i)
#         otusam1 = otuwsa_od.get(i)[0]
#         otusam2 = otuwsa_od.get(i)[1]
#         if fp <= 0.05:
#             qf+=1

#         if tp <= 0.05:
#             qt+=1

#     print("fisher, chi ={}".format(qf),"ttest ={}".format(qt),  sep="\t", end="\n", file=log )


idfs.close()
    # hmn = [otu_table.get_value_by_ids('denovo228',i) for i in u]
    # i = otu_table.iter_data()
        # ind = sumotudict.keys().index(i)
        # indm = ind + 1 in enumerate(zip(otu_list[:-1], otu_list[1:]))
        # sum2 = sumotudict.values()[indm]
        