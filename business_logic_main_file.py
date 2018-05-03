#!/usr/bin/env python3

"""
Program:        chromosome15_gene_browser
File:           business_logic_main_file.py
Version:        V1.0
Date:           17.04.2018
Function:       Provide a simple file for a front end to work with. See business_logic.py file
                for break down of functions
Copyright:                  Sarah Griffiths
Author:                     Sarah Griffiths
project_collaborators:      Archana Patil, Sergio Ramalli, Fabio Biond
Address:                    Bioinformatics, BBK, London
-------------------------------------------------------------
This program is released under the General use (GNU) Public License
-------------------------------------------------------------
Description:
==============
This programme accepts requests from a web page to query the data base Chromosome15
and processes and returns requested data in Python Variables - see user guide for description
of variables
---------------------------------------------------------------
Usage:
===========
Chromosome15_gene_browser
-----------------------------------------------------------------
Revision History:
=================
V1.0 17.04.2018     Original    BY: Sarah Griffiths
v2.0 28.04.2018     Original    BY: Sarah Griffiths
"""
#import business_logic
#import archanapymysql as DB

def getAllEntryData(key):

    XYZ to retrieve the data. 

    parsedSequence = BL.ParseSequence(new)
    codingRegion = BL.codingRegion(start,end,parsedSequence)
    mrnaSequence = BL.translate(codingRegion)
    splitSequence = BL.CodonSequence(mrnaSequence)
    translatedAndAligned = BL.alignseq(splitSequence)
    justAminoAcids = BL.translatedSequence(splitSequence)#check this
    codonFrequency = BL.codonFreq(splitSequence)# need to edit to incorporate total frequencies
    restrictionEnzymeCutSites = BL.restrictionEnzyme('ttgtc', start, end, parsedSequence) #returned as dictionary

