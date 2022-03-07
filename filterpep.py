#!/bin/env python
from sys import argv
# coding: utf-8

infile = open(argv[1])
for line in infile :
    line2 = line.strip().split("\t")
    ID = line2[0]
    seq = line2[1]
    Ccount = seq.count("C")
    Scount = seq.count(".")
    if Scount ==0 and Ccount <=5 :
        print(ID,seq,"pass",sep="\t")
    elif Scount == 0 and  Ccount >= 6 :
        print(ID,seq,"Cys-Rich",sep="\t")
    else :
        continue

