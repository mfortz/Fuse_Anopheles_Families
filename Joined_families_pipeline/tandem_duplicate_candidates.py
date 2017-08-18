#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

tandem_duplicate_candidates.py

Created on Wed Jul 26 19:14:20 2017

@author: mfortz
"""

from __future__ import division

import sys        
from deriveGeneStructure import deriveGeneStructure

all_gene_file = open(sys.argv[1],'r').readlines()
data, geneOrder = deriveGeneStructure(all_gene_file)

min_species_spanned_by_f1 = int(sys.argv[2])
if min_species_spanned_by_f1 < 2:
    print("Error. Max size must be >= 2.")
    sys.exit()
    
ratioIn = float(sys.argv[3])

#output file
candidatesOut = open(sys.argv[4], 'w+')
candidatesOut.write("F1 \t F2_gene \n")

nonEmptySpecies = {}
#f1set = []

    
for fam in data.familiesDict:
    nonEmptySpecies[fam] = 0    
    currentFamily = data.familiesDict[fam]
    
    for species in data.speciesDict:
        if len(currentFamily.familyMembersDict[species]) > 0:
            nonEmptySpecies[fam] += 1



Result = {}

for gene in data.genesDict:
    currentGene = data.genesDict[gene]
    currentFamName = currentGene.family
    
    #look for f1 candidates, ie less than max size
    if nonEmptySpecies[currentFamName] >= min_species_spanned_by_f1:
        
        #load an ordered list of genes on the current contig and species
        currentFragment = geneOrder[currentGene.species][currentGene.ctg]   
        positionOnCtg = currentFragment.index(currentGene) 
        
       
        #adds the neighbour to the left or to the right
        #offset = 0 if left nbr, 1 if right nbr
        def addNbr(offset):  
            
            #takes care of errors arising from the gene being at the end of its contig
            if offset == 0:
                if positionOnCtg == 0:
                    return
                else:
                    nbr = currentFragment[positionOnCtg - 1]
            elif offset == 1:
                try:
                    nbr = currentFragment[positionOnCtg + 1]
                except IndexError:
                    return
    
            #aif f1 != f2 and size_f1 <= size_f2, add to the appropriate list
            if currentFamName != nbr.family and nonEmptySpecies[currentFamName] <= nonEmptySpecies[nbr.family]:
                try:
                    Result[(currentFamName,nbr.family)].append((currentGene.gene,nbr.gene))
                except KeyError:
                    Result[(currentFamName,nbr.family)] = []
                    Result[(currentFamName,nbr.family)].append((currentGene.gene,nbr.gene))
        
        
        #add both left and right neighbours to Result
        addNbr(0)
        addNbr(1)


        
"""Filter result for pairs such that the ratio of syntenic order to max_size is above threshold"""
doneMap = {}

for f1,f2 in Result:
    
    """Filter for repeated entries, ie instances of (f1,f2) and (f2,f1)"""
    try:
        doneMap[(f2,f1)]
        continue
    except KeyError:
        pass
    
    observedSynteny = len(Result[(f1,f2)])
    f1size = nonEmptySpecies[f1]
    ratio = observedSynteny/f1size
    
    if ratio >= ratioIn:
        doneMap[(f1,f2)] = None
        doneMap[(f2,f1)] = None
        
        candidatesOut.writelines([f1,'\t',f2,'\n'])
        
        

        
        
        



        
        
        
        
        
        
    
    