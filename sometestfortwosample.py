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
from functools import wraps

inp = argv[1]
ID = argv[2]

class Dict():



    def __init__(self):
        self.idfs = open(ID, "r")
        self.otu_table = biom.load_table(inp)


    def sample_list(self):
        a = self.ftable()
        sample_list = a.ids(axis='sample') #all sample ids
        return sample_list

    def otu_list(self):
        a = self.ftable()
        otu_list = a.ids(axis='observation') #all obs ids(otus)
        return otu_list

    def ftable(self):
        '''
        Filter biom table by sample Id from mapping txt file.
        '''
        al = self.idfs.read()
        al_l =[line.rstrip() for line in al.split("\t")]
        al_ll = []
        for a in al_l:
            if re.search(r'\D\Z',a):
                sample_list = self.otu_table.ids(axis='sample')
                for i in sample_list: 
                    if i.startswith(a):
                        al_ll.append(i)
            else:
                al_ll.append(a)
        
        ftable = self.otu_table.filter(al_ll,axis='sample', inplace=False)
        return ftable

    def iddict(self):
        '''
        Ordered dict with groupped otu sample: key - master, value - master replication
        '''
        samplesame_list = []
        sample_list = self.sample_list()
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

    def delete_minors(self): 
        otu_del = [ d[0] for d in self.summoldict.items() if self.summoldict.get(d[0])[0] >= 6 and summoldict.get(d[0])[1] ]
        delete_minors = ftable.filter(otu_del,axis='observation', inplace=False)   
        return delete_minors

    def sumotudict(self):
        '''
        Create OrderedDict odject
        key = OTU Id
        value = list of sum for all repeats in samle order
        '''
        sumotudict = OrderedDict()
        for q in self.otu_list():
            suml = []
            a = self.iddict()
            for i in a.values():
                ftable = self.ftable()
                w = [ftable.get_value_by_ids(q, a) for a in i]
                sum_otu = sum(w)
                suml.append(sum_otu)
            sumotudict.update({q:suml})
        return sumotudict

    def sumalldict(self):
        '''
        Function create OrderedDict object.
        key = sample Id
        value = sum of all otu for this sample
        '''
        sumalldict = OrderedDict()
        a = self.sample_list()
        leng = len(a)
        lenglist = range(0, leng)
        a = self.iddict()
        keyids = a.keys()
        zipl = zip(keyids, lenglist)
        for w in zipl:
            a = w[1]
            s = self.sumotudict()
            r = s.values()
            y =[t[a] for t in r]
            sum_y = sum(y)
            sumalldict.update({w[0]:sum_y})
        return sumalldict

    def otsa_od(self):
        '''
        Ordered dict with a pair lists in values for tested samples.
        '''
        otsa_od = OrderedDict()
        for i in self.otu_list:
            dpartbef = [a for a in self.iddict.values()]
            dpart1 = [self.ftable.get_value_by_ids(i, a) for a in dpartbef[0]]
            dpart2 = [self.ftable.get_value_by_ids(i, a) for a in dpartbef[1]]
            otsa_od.update({i:[dpart1,dpart2]})
        return otsa_od

    def sometests(self, pdict):

        '''
        This function consists of statistical tests, that calculate p-value for our data.
        This is my  shiny shit castle of crap, which really lacks some order.
        Return orderdicts with pairs like  otu:p-value: 
        fisher/chi2(if sum from some sample more than 5)[0]
        ttest [1]
        kruskal[2]
        '''

        otherdict = OrderedDict()
        a = self.sumalldict()
        summall = a.values()
        for i in self.sumotudict.items():
            other = map(operator.sub, summall , i[1])
            otherdict.update({i[0]:other})

        leno = len(otu_list)
        fj=0
        pdict = OrderedDict()
        for i in self.sumotudict.items():
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
            a = self.otuwsa_od.get(i)[0]
            b = self.otuwsa_od.get(i)[1]
            p = ttest(a,b)[1]
            if p != p:
                p = 1
            tdict.update({i:p})

            sys.stdout.write('\r')
            sys.stdout.write("ttest {}/{}".format(j,leno))
            sys.stdout.flush()

        mdict = OrderedDict()
        mj=0
        for i in otu_list:
            mj+=1
            a = self.otuwsa_od.get(i)[0]
            b = self.otuwsa_od.get(i)[1]
            if np.sum(a) == 0 or np.sum(b) == 0:
                p = "nan"
            else:
                p = man(a,b)[1]
            mdict.update({i:p})
            sys.stdout.write('\r')
            sys.stdout.write("man-y {}/{}".format(mj,leno))
            sys.stdout.flush()

        kdict = OrderedDict()
        kw=0
        for i in otu_list:
            kw+=1
            a = self.otuwsa_od.get(i)[0]
            b = self.otuwsa_od.get(i)[1]
            p = wilc(a,b)[1]
            kdict.update({i:p})
            sys.stdout.write('\r')
            sys.stdout.write("kruskal {}/{}".format(kw,leno))
            sys.stdout.flush()
        
        sys.stdout.write('\r')
        
        return pdict, tdict, kdict, mdict

    def stat_filter(self):
        fdd = [d[0] for d in self.pdict.items() if d[1] <= 0.05]
        stat_filter = self.ftable.filter(fdd, axis='observation', inplace=False)
        return stat_filter

    def fbiom_differ(self):
        '''division and subtration parts with each other'''
        fddiv = OrderedDict()
        part = OrderedDict()
        a = self.sumotudict()
        for d in self.sumotudict.items():
            key = d[0]
            summ_i = self.sumalldict.values()
            summ_1 = summ_i[0]
            summ_2 = summ_i[1]
            val_w = d[1]
            pval_w_1 = val_w[0]/summ_1
            pval_w_2 = val_w[1]/summ_2
            part.update({key:[pval_w_1, pval_w_2]})
        for d in part.items():
            key = d[0]
            w = d[1]
            if w[0] == 0:
                div = 0.0000000001 / w[1]
            elif w[1] == 0:
                div = w[0] / 0.0000000001
            else:
                div = w[0] / w[1]
            fddiv.update({key:div})
        fdsub = OrderedDict()
        for d in part.items():
            key = d[0]
            w = d[1]
            sub = w[0] - w[1]
            fdsub.update({key:sub})
        return fddiv, fdsub, part
    
 

class Flow(Dict):

    def __init__(self):
        Dict.__init__(self)

    def main_pseudo(self):
        '''hop, hop pseudozeros and minor otu eliminating'''
        table_del = delete_minors(self.sumotudict, self.ftable)
        print(table_del)
        return table_del

    def main_old(self):
        pdict, tdict, kdict, mdict = self.sometests(self)
        fd = self.fbiom_differ(pdict, self.ftable) 
        return fd

    def deletions(self):
        '''start doing monkey job again for new filtered OTU table'''

        sys.stdout.write('\r')
        sys.stdout.write("Doing divisions and subtrations for greater good.")
        sys.stdout.flush()


        fd = self.main_old()
        fotu_table = fd.ids(axis='sample')()
        otu_list = fd.ids(axis='observation')()
        s_iddict = self.iddict(fotu_table) 
        otuwsa_od = self.otsa_od(s_otu_list,s_iddict) 
        f_sumotudict = self.summotu_od(s_otu_list,s_iddict) 
        fddiv, fdsub, part = self.fbiom_differ(sumotudict, f_sumotudict, sumalldict)
        return otu_list, fddiv, fdsub, part



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

# with open("log.txt", "a") as log:
#     for i in otu_list:
#         pdict = pvalues[0]

#         fp = pdict.get(i)
#         tp = tdict.get(i)
#         otusam1 = otuwsa_od.get(i)[0]
#         otusam2 = otuwsa_od.get(i)[1]
#         if fp <= 0.05:pe
#             qf+=1


#         if tp <= 0.05:
#             qt+=1




sys.stdout.write('\r')
sys.stdout.write("Prepare biom table for analysis")
sys.stdout.flush()


'''mainflow'''

a = Flow()
# print(a.sample_list())
print(a.otu_list())

# with open("subdel.txt", "w") as out:
#     print("otu", "division", "subtration","part_1","part_2", sep="\t", end="\n", file=out)
#     a = Flow()
#     otu_list, fddiv, fdsub, part = a.deletions()
#     for i in otu_list:
#         p1 = fddiv.get(i)[0]
#         p2 = fddiv.get(i)[1]
#         print(i, div,sub, p1, p2, sep="\t", end="\n", file=out)

#     print("fisher, chi ={}".format(qf),"ttest ={}".format(qt),  sep="\t", end="\n", file=log )


    # hmn = [otu_table.get_value_by_ids('denovo228',i) for i in u]
    # i = otu_table.iter_data()
        # ind = sumotudict.keys().index(i)
        # indm = ind + 1 in enumerate(zip(otu_list[:-1], otu_list[1:]))
        # sum2 = sumotudict.values()[indm]
        

''' Some fancy plots. Try to realize this part on pandas libraries '''

