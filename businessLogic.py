
#!/usr/bin/env python3

"""
Program:        chromosome15_gene_browser
File:           buseinss_layer.py
Version:        V2.1
Date:           03.05.2018
Function:       Obtain data from Chromosome15 DB, process and return output in JSON format.
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
and processes and returns requested data in Python Variables
---------------------------------------------------------------
Usage:
===========
Chromosome15_gene_browser
-----------------------------------------------------------------
Revision History:
=================
V1.0 24.03.2018     Original    BY: Sarah Griffiths
V2.0 17.04.2018     Re-Work 	BY:Sarah Griffiths
v2.1 03.05.2018	    Fixed a bug	BY: Sarah Griffiths
v2.2 11.05.2018     added cds join codon BY:Sarah Griffiths
"""
#***********************************************************
#import libraries
#from urllib import request
import sys
import re
#import archanapymsql

#************************************************************
""" Classes and functions"""
#could make this and the one below a class.
#could store codons in seperate file and import module codons?

class BasicDBExtraction(object):
      """ Class containing all functions required to extract basic info from DB with no computations the returned values can then be used in calculation - you must call these function and classes first
"""

      def BasicInfo(gene_name):
            """input = Gene_name
                  Output = accession number, NCBI identifier,
                  chromosome location, protein product name -
                  dictionary containing these items
				  17.04.2018 		By:SG
            """
                  
            BasicInfoDictionary = Gene (gene_name)
            
            locals().update(BasicInfoDictionary)

            return (accessionNumber)
            return(NCBIIdentifier)
            return (chromosomeLocation)
            return (proteinProductName)

      def SequenceExtraction(gene_name):
            """input = Gene_name
                  Output = CDS start position, CDS end
                  position, gDNASequence and CDS -
                  dictionary containing these items
				  17.04.2018 		By:SG
            """
                  
            SequenceDictionary = Sequence(gene_name)
            
            locals().update(SequenceDictionary)

            

            return (CDSStartposition)
            return(CDSendposition)
            return (gDNASequence)
            return (CDS)

      def RestrictionEnzyme(RE_name):
            """input = RE_name
                  Output = RE sequence
				  17.04.2018 		By:SG
            """
                  
            REDictionary = Restriction_enzyme(RE_name)
            
            locals().update(REDictionary)

            
      
            return (REsequence)
    
                        


def ParseSequence(origin):
    """ Parse the origin sequence returned from chromomse15_DB
    Input:  origin     --- string returned by pymysql in back end 'origin' is name of sequence in genbank file
    Return: sequence   --- parsed origin so now this is just nucleic acid sequence
    24.03.2018         By: SG
    17.04.2018		By:SG
    """
    count = 0
    sequence = ''
    for item in origin:
        if count == 0:
            count += 1
        elif count % 7 == 0:
            count += 1
        else:
            count += 1
            sequence += item
    return (sequence)


def codingRegion(start,end,sequence):
    """Return the coding region based on location provided # note to self may need to parse the location to get start&end
    Input : (start, end) --- integers returned by pymysql in back end on request from front end
    Return: gene         --- string of nucleic acids that code for gene
    24.03.2018         By: SG
    """
    count = 0
    gene = ''
    start -= 1
    for nucleotide in sequence:
      if count >= start and count <=end:
        gene += nucleotide
        count += 1
      else:
        count += 1
    return (gene)


def translate(code):
    """
    translate dna seq to RNA sequence
    Input: sequence         --- nucleic acid sequence
    Return: mrna            --- t's replaced with u's
    24.03.2018  By: SG
    """
    mrna = code.replace('t', 'u')
    return mrna

def CodonSequence(rna):
    """
    translate RNA sequence in to Amino acid sequence
    Input:  rna  --- mrna string
    Return: codon_sequence --- nucleic acids split in to 3's for amino acid sequence
    24.03.2018 By: SG
    """


    length = 3
    codon_sequence = [rna[i:i + length] for i in range(0, len(rna), length)]

    return (codon_sequence)
# for a test
def translatedSequence(codon_sequence):
      """ A test to check translation of sequence is the same as CDS in records #ammend once DB layer done
      Input: codon_sequence
      Return: Amino acid sequence as a stringand 'Match or not match'"""

      codons = {"uuu": "F", "uuc": "F", "uua": "L", "uug": "L",
              "ucu": "S", "ucc": "S", "uca": "S", "ucg": "S",
              "uau": "Y", "uac": "Y", "uaa": "STOP", "uag": "STOP",
              "ugu": "C", "ugc": "C", "uga": "STOP", "ugg": "W",
              "cuu": "L", "cuc": "L", "cua": "L", "cug": "L",
              "ccu": "P", "ccc": "P", "cca": "P", "ccg": "P",
              "cau": "H", "cac": "H", "caa": "Q", "cag": "Q",
              "cgu": "R", "cgc": "R", "cga": "R", "cgg": "R",
              "auu": "I", "auc": "I", "aua": "I", "aug": "M",
              "acu": "T", "acc": "T", "aca": "T", "acg": "T",
              "aau": "N", "aac": "N", "aaa": "K", "aag": "K",
              "agu": "S", "agc": "S", "aga": "R", "agg": "R",
              "guu": "V", "guc": "V", "gua": "V", "gug": "V",
              "gcu": "A", "gcc": "A", "gca": "A", "gcg": "A",
              "gau": "D", "gac": "D", "gaa": "E", "gag": "E",
              "ggu": "G", "ggc": "G", "gga": "G", "ggg": "G", }
      aminoacidSequence = ''

      for i in (codon_sequence):
            if i in codons:
                  aminoacidSequence += codons[i]
            else:
                  aminoacidSequence += '?'
      return(aminoacidSequence)            
            

      

      
      

def alignseq(aminoacid):

    """
    Align nucleic acid and amino acid sequence
    Input: amino acid sequence --- either from DB or calculated seq
    Return: Aligned sequence
    24.03.2018  By: SG
    """
    codons = {"uuu": "F", "uuc": "F", "uua": "L", "uug": "L",
              "ucu": "S", "ucc": "s", "uca": "S", "ucg": "S",
              "uau": "Y", "uac": "Y", "uaa": "STOP", "uag": "STOP",
              "ugu": "C", "ugc": "C", "uga": "STOP", "ugg": "W",
              "cuu": "L", "cuc": "L", "cua": "L", "cug": "L",
              "ccu": "P", "ccc": "P", "cca": "P", "ccg": "P",
              "cau": "H", "cac": "H", "caa": "Q", "cag": "Q",
              "cgu": "R", "cgc": "R", "cga": "R", "cgg": "R",
              "auu": "I", "auc": "I", "aua": "I", "aug": "M",
              "acu": "T", "acc": "T", "aca": "T", "acg": "T",
              "aau": "N", "aac": "N", "aaa": "K", "aag": "K",
              "agu": "S", "agc": "S", "aga": "R", "agg": "R",
              "guu": "V", "guc": "V", "gua": "V", "gug": "V",
              "gcu": "A", "gcc": "A", "gca": "A", "gcg": "A",
              "gau": "D", "gac": "D", "gaa": "E", "gag": "E",
              "ggu": "G", "ggc": "G", "gga": "G", "ggg": "G", }
    upper = ''
    lower = ''
    for acid in aminoacid:
        if acid in codons:
            
            code = codons[acid]
            upper += acid

            lower += '--'
            lower += code
        else:
            pass
    return (upper, lower)





def codonFreq(nuc_sequence):
    """ Calculate the codon frequency and returning a dictionary- @need to edit this so it returns a frequency once we have total
     Input: codons from CDS of particular gene  --- returned from coding region
     Returns: Codon frequency in a dictionary
     24.03.2018  By: SG"""
    codons = {"uuu": "F", "uuc": "F", "uua": "L", "uug": "L",
              "ucu": "S", "ucc": "s", "uca": "S", "ucg": "S",
              "uau": "Y", "uac": "Y", "uaa": "STOP", "uag": "STOP",
              "ugu": "C", "ugc": "C", "uga": "STOP", "ugg": "W",
              "cuu": "L", "cuc": "L", "cua": "L", "cug": "L",
              "ccu": "P", "ccc": "P", "cca": "P", "ccg": "P",
              "cau": "H", "cac": "H", "caa": "Q", "cag": "Q",
              "cgu": "R", "cgc": "R", "cga": "R", "cgg": "R",
              "auu": "I", "auc": "I", "aua": "I", "aug": "M",
              "acu": "T", "acc": "T", "aca": "T", "acg": "T",
              "aau": "N", "aac": "N", "aaa": "K", "aag": "K",
              "agu": "S", "agc": "S", "aga": "R", "agg": "R",
              "guu": "V", "guc": "V", "gua": "V", "gug": "V",
              "gcu": "A", "gcc": "A", "gca": "A", "gcg": "A",
              "gau": "D", "gac": "D", "gaa": "E", "gag": "E",
              "ggu": "G", "ggc": "G", "gga": "G", "ggg": "G", }

    codonfrequency = {}
    for codon in nuc_sequence:
        if codon in codonfrequency:
                codonfrequency[codon] += 1
        else:
                codonfrequency[codon] = 1
        for i in codons:
            if i not in codonfrequency:
                codonfrequency[i] = 0
    return (codonfrequency) 


def totalCodonFreq(seq):
    """ Calculate codon frequency across all coding regions @ need to think of a way to deal with outliar codons?
    Input: all coding regions                   ---
    Returns:  Codon frequency in a dictionary --- back end will then need to make this in to a table to store data for comparison later
    24.03.2018  By: SG"""

    #split in to 3's and then calculate


    n = 3

    all_mrna = seq.replace('t', 'u')
    all_mrna = [all_mrna[i:i + n] for i in range(0, len(all_mrna), n)]

    total_freq = {}
    for codon in all_mrna:
        if codon in total_freq:
            total_freq[codon] += 1
        else:
            total_freq[codon] = 1
    for i in codons:
        if i not in freq:
            total_freq[i] = 0
    return (total_freq)


def restrictionEnzyme(shortSequence, cds_start, cds_end,sequence):
    """ Identifies restriction enzyme cut sites by looking for palindromes in sequence and signals if they are inside or outside coding regions
    Input: sequence --- restriction enzyme site selected from DB list or input yourself
           cds_start --- start of CDS which can be retrieved from DB @SG will need to retrieve these with pymysql once back end complete - for now just use dummys
           cds_end  --- end of CDS which can also be retrieved from the DB
    Output: dictionary containing start and end locations of entire palindrome and whether or not they ae in the coding region
     # note to self - perhaps give variable names to match.start() & end """
    reverse = shortSequence[::-1]
    pattern = re.compile(r'(' + shortSequence + ').+?(' + reverse + ')')
    sites = {}
    for match in pattern.finditer(sequence):
        if match.start() > cds_start and match.end() < cds_end:
            sites[match] = (match.start(), match.end(), 'in coding region')

        else:
            sites[match] = (match.start(), match.end(), 'not in coding region')

    for i in sites:
        print(i, sites[i])
    return(sites)   

def cdsJoin(join):
      """ function to deal with format of CDS data
      Input: Join from database 
      Output: start and end variables
      """
      
      p = re.compile(r'<(\d+).+>(\d+)')
      match = p.search(join)

      start = match.group(1)
      end = match.group(2)
      return (start,end)

