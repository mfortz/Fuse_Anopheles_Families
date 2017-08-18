#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plotHelper.py

Methods to plot silhouette distributions. Also prints mean and standard 
deviation of the distribution. Usage:

	$ python plotHelper.py <module> <input>

Modules:
	1. 'all_families' - silhouette distribution on all families and all genes.
		Input: all silhouette file.
	2. 'joined' - change in silhouette distribution after joining family pairs.
		Input: joined families silhouette.

Created on Thu Jul 27 20:34:06 2017

@author: mfortz
"""

import sys
import matplotlib, matplotlib.pyplot as plt
import numpy as np

#inputs
mode = sys.argv[1]
fileIn = open(sys.argv[2],'r').readlines()



"""Plot distribution of silhouette difference after joining families
Input: joined_silhouette"""
def plotJoinedFamilies(fileIn):
    differenceList = []
    
    #list scores as float
    for l in fileIn[1:]:
        x = l.split()
        d = float(x[7])
        differenceList.append(d)
    
    #statistics
    mean = np.mean(differenceList)
    stddev = np.std(differenceList)
    shortmean = "%.4f" % mean
    shortstd = "%.4f" % stddev
    print("Mean =", shortmean,". Standard deviation =",shortstd)
        
    #plot
    plt.rcdefaults()
    fig, ax = plt.subplots()

    binArray = np.arange(-1.0, 1.05, 0.1)
    
    fig.set_size_inches(6,6)
    ax.set_xlim([-1,1])
    ax.set_yscale('log', nonposy='clip')
    
    ax.hist(differenceList, bins=binArray, color='#207de3', edgecolor='none')
    ax.set_xlabel('Silhouette change')
    ax.set_ylabel('Number of family pairs')
    ax.set_title('Distibution of silhouette difference')
    
    plt.show()



""" Plot distribution of silhouette on families and genes
Input: all_fam_silhouette"""    
def plotGeneFamSilhouette(fileIn):
    geneSilhouette = []
    famSilhouette = []
    
    #list scores as loat
    for l in fileIn[1:]:
        x = l.split()
        
        #add ave family silhouette
        famScore = float(x[5]) 
        famSilhouette.append(famScore)
        
        #add gene silhouettes
        i = 7
        while i < len(x):
            geneScore = float(x[i])
            geneSilhouette.append(geneScore)
            i += 2
     
    #statistics
    meanf = np.mean(famSilhouette)
    stddevf = np.std(famSilhouette)
    shortmeanf = "%.4f" % meanf
    shortstdf = "%.4f" % stddevf
    print("Families: \t Mean =", shortmeanf,". Standard deviation =",shortstdf)

    meang = np.mean(geneSilhouette)
    stddevg = np.std(geneSilhouette)
    shortmeang = "%.4f" % meang
    shortstdg = "%.4f" % stddevg
    print("Genes: \t Mean =", shortmeang,". Standard deviation =",shortstdg)

    #plots        
    plt.rcdefaults()
    fig, (ax1,ax2) = plt.subplots(1,2)
    
    fig.set_size_inches(5,20)

    binArray = np.arange(-1.0, 1.05, 0.1)
    
    ax1.set_xlim([-1,1])
    #ax1.set_ylim([0,10**6])
    ax1.set_yscale('log', nonposy='clip')
    
    ax2.set_xlim([-1,1])
    #ax2.set_ylim([0,10**7])
    ax2.set_yscale('log', nonposy='clip')

    ax1.hist(geneSilhouette, bins=binArray, color='#1165c1', edgecolor='none')
    ax1.set_xlabel('Silhouette')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Silhouette of genes')
    
    ax2.hist(famSilhouette, bins=binArray, color='#1165c1', edgecolor='none')
    ax2.set_xlabel('Silhouette')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Silhouette of families')
    
    fig.subplots_adjust(wspace=.5)
    plt.show()        
        
        
"""Run modules"""
if mode == 'all_families':
    plotGeneFamSilhouette(fileIn)
elif mode == 'joined':
    plotJoinedFamilies(fileIn)
    
    
    
    