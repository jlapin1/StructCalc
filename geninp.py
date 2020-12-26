#! /usr/bin/python3

import numpy as np
import os
import sys

def tensor(Iso):
    # Paramasivan et al. 2018
    sig11 = 0.96*Iso + 99.7
    sig22 = 1.16*Iso - 55.1
    sig33 = 0.88*Iso - 44.5
    return [sig33, sig22, sig11]

def readft(filename):
    f = open(filename, "r")
    cnyq=0;rnyq=0;cnyqd=1;rnyqd=1;refc=0;refr=0
    for line in f:
        ln = line.split()
        if len(ln)>0:
            if ln[0]=='-xSW':
                cnyq = float(ln[1])
                rnyq = float(ln[3])
            if ln[0]=='-xOBS':
                cnyqd = float(ln[1])
                rnyqd = float(ln[3])
            if ln[0]=='-xCAR':
                refc = float(ln[1])
                refr = float(ln[3])
    return cnyq, cnyqd, rnyq, rnyqd, refc, refr

fread = {
    'inp':False,
    'seq':False,
    'iso':False,
    'ft':False
}

# Set working directory
dir = input("Type path to project directory\n>>> ")
os.chdir(dir)

# Read in information from all relevant input files
ifiles = os.listdir("raw/")
for m in ifiles:
    if m[-3:]=='inp':
        with open("raw/"+m,"r") as f:
            print("\nReading raw targets in from %s"%(m))
            f.readline()
            sz = int(f.readline())
            rtargs = np.zeros((sz, 3))
            for n in range(sz):
                rtargs[n] = f.readline().split(",")
        fread['inp'] = True
    if m[-3:]=='seq':
        with open("raw/"+m,"r") as f:
            print("\nReading sequence in from %s"%(m))
            f.readline()
            resraw = f.readline().strip()
        fread['seq'] = True
    if m[-3:]=='iso':
        with open("raw/"+m,"r") as f:
            print("\nReading isotropic shifts in from %s"%(m))
            f.readline()
            sz2 = int(f.readline())
            iso = np.zeros((sz2))
            for n in range(sz2):
                iso[n] = float(f.readline())
        fread['iso'] = True
    if m[-3:]=='com':
        print("\nReading Nyquist, carrier, and ref. frequencies from %s"%(m))
        xnyq, xcar, ynyq, ycar, xref, yref = readft("raw/"+m)
        print("Values found\nX-Nyquist: %.1f\nY-Nyquist: %.1f\nX-Carrier: %.3f\n\
Y-Carrier: %.3f\nX-Ref.: %.1f\nY-Ref.: %.1f"%(xnyq, ynyq, xcar, ycar, xref, yref))
        fread['ft'] = True
        ans = str(input("\nAre you OK with parameters read in from .ft file? (y/n)\n>>> "))
        if ans=='n':
            ans = input("Input comma separated values Xnyq,Ynyq,Xcar,Ycar,Xref,Yref\n>>> ")
            [xnyq, ynyq, xcar, ycar, xref, yref] = [float(m) for m in ans.split(",")]

# Quality control
if fread['inp']==False:
    sys.exit("\nError: No .inp file found")

# ft.com
if fread['ft']==False:
    print("\nNo .com file found")
    ans = str(input("Input comma separated values Xnyq,Ynyq,Xcar,Ycar,Xref,Yref\n>>> "))
    [xnyq, ynyq, xcar, ycar, xref, yref] = [float(m) for m in ans.split(",")]

# Create isotropic shifts if no file was found
if fread['iso']==True:
    sigmas = np.array([tensor(m) for m in iso])
else:
    sigma = [64, 77, 222]
    sigmagly = [41, 64, 215]
    print("\nNo isotropic shifts file found\nDefaulting to\n\t[%d %d %d] for glycine\n\t[%d %d %d] for all else"%(sigmagly[0],sigmagly[1],sigmagly[2],sigma[0],sigma[1],sigma[2]))
    ans = input("\n\nAre you ok with these tensor parameters (y,n)?\n>>> ")
    if ans.upper() != 'Y':
        sigmagly = [float(n) for n in input("\n\nInput 3 space separated values for glycine (# # #)\n>>> ").split()]
        sigma = [float(n) for n in input("\nInput 3 space separated values for all other residues (# # #)\n>>> ").split()]
    iso = [];avg = []
    sigmas = np.zeros(( sz, 3 ))
    for i,j in enumerate(resraw):
        three = sigmagly if j=='G' else sigma
        iso.append(np.mean(three))
        for k,l in enumerate(three):
            sigmas[i,k] = l*(iso[-1]/np.mean(three))

S0 = float(input("\nInput order parameter So (~0.8)\n>>> "))

length = rtargs.shape[0]
CS = np.zeros((length))
NH = np.zeros((length))
CH = np.zeros((length))
for m in range(length):
    CS[m] = xcar*( ((rtargs[m,0]-iso[m]))/(-S0/2) + iso[m] )
    NH[m] = rtargs[m,1]/(S0/2)
    CH[m] = rtargs[m,2]/(S0/2)

with open("input/realtarg.csv","w") as f:
    f.write("# Header\n")
    f.write("%d\n"%(length))
    for m in range(length):
        f.write("%c,%.2f,%.2f,%.2f\n"%(resraw[m], CS[m], NH[m], CH[m]))
print("\nWrote realtarg.csv to direcotry input/")

with open("input/sigin.csv","w") as f:
    f.write("# Header\n")
    for m in sigmas:
        f.write("%.2f,%.2f,%.2f\n"%(m[0],m[1],m[2]))
print("Wrote sigin.csv to directory input/")

with open("output/angles.csv","w") as f:
    f.write("So,%.3f\n"%(S0))
print("Wrote order parameter So to output/angles.csv")
