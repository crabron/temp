#! /usr/bin/env python

from sys import argv
import numpy as np
import biom
import pandas
import h5py


inp = argv[1]


def main(inp):
    otu_table = biom.load_table(inp)
    np.float64 == np.dtype(float).type
    a = otu_table.to_dataframe()
    print(a)



if __name__ == "__main__":
    main(inp)
