#! /usr/bin/env python

"""
Work only for couples. Ignor groups. Relise of preveous code, but more kinda introvertic.
"""

from __future__ import division, print_function
from scipy.stats import fisher_exact as fisher
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


inp = "random_first.biom"
idfs = open("IDtestlist.txt", "r")
otu_table = biom.load_table(inp)
ftable = filter_otu(otu_table, idfs)
sample_list = ftable.ids(axis='sample')
otu_list = ftable.ids(axis='observation')

gs = group_sample(sample_list)

key_id = [re.split(r'.\d',i[0])[0] for i in gs]
zipl = zip(key_id,gs)
iddict = OrderedDict(zipl)

sumotudict = summotu_od(otu_list,iddict)
sumalldict = summall_od(otu_list,iddict, sumotudict)


idfs.close()

otherdict = OrderedDict()
summall = sumalldict.values()
for i in sumotudict.items():
    other = map(operator.sub, summall , i[1])
    otherdict.update({i[0]:other})

# bonff = 
pfishdict = OrderedDict()
for i in sumotudict.items():
    sum = i[1]
    other = otherdict.get(i[0])
    table = np.array([sum,
                    other])
    # if sum[1] <= 5 or sum[0] <= 5 
    p = fisher(table)[1]
    pfishdict.update({i[0]:p})
    print(p)




idfs.close()
# hmn = [otu_table.get_value_by_ids('denovo228',i) for i in u]
# i = otu_table.iter_data()
    # ind = sumotudict.keys().index(i)
    # indm = ind + 1 in enumerate(zip(otu_list[:-1], otu_list[1:]))
    # sum2 = sumotudict.values()[indm]

