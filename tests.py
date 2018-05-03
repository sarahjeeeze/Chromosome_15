#!/usr/bin/env python3

"""
Program:        chromosome15_gene_browser
File:           tests.py
Version:        V1.0
Date:           28.04.2018
Function:       Tests to identify problems with business_logic or data errors
                
Copyright:                  Sarah Griffiths
Author:                     Sarah Griffiths
project_collaborators:      Archana Patil, Sergio Ramalli, Fabio Biond
Address:                    Bioinformatics, BBK, London
-------------------------------------------------------------
This program is released under the General use (GNU) Public License
-------------------------------------------------------------
Description:
==============
Tests business logic layer
---------------------------------------------------------------
Usage:
===========
Chromosome15_gene_browser
-----------------------------------------------------------------
Revision History:
=================
v2.0 28.04.2018     Original    BY: Sarah Griffiths

"""
def testTranslation(translated):
    """ run a test to see if translated sequence contains expected number of AA
        or if it had lost any due anomolous codons. 
    Input: translated  -- translated sequence
            start -- cds start
            end -- cds end
     Return:  Success/more than/ less than     
    """
    expected = (end - start) / 3
    count = 0
    for aminoAcid in translated:
        count += 1
    if count == expected:
        print("Success")
    elif count < expected:
        print("less than expected number of AA")
    else:
        print("more than expected number of AA")

    
        
def compareTranslationAndDatabaseSequence(justAminoAcids, expected):
    """ compares translation with expected translation either pasted in or retrieved from DB
        Input : justAminoAcids -- variable returned from running getAllEnteyData
                expected -- either retrieved from DB or pasted in
        Return: "match" or "no match"
    """

    if justAminoAcids == expected:
        print('match')
    else:
        print('no match')
        
def expectedCodonFrequencies(codonFrequency):
    """ counts number of dictionary items and alerts if greater than 64
    Input: codonFrequency - output from running getAllEntryData
    Return: "expected number of codons" or "more than expected - must include rare codons"
    """
    count = 0
    for key in codonFrequency:
        count += 1
    if count >= 64:
        print('more than expected - must include rare codons')
    else:
        print("expected number of codons")
        

                
            
