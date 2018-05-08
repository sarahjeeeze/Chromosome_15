import businessLogic as BL

#dummy data
with open ('seq.txt','r') as f:
    sequence = f.read().split()

print(sequence)
import experiment

#dummy 
start = 3360
end = 4300

parsedSequence = BL.ParseSequence(sequence)
codingRegion = BL.codingRegion(start,end,parsedSequence)
mrnaSequence = BL.translate(codingRegion)
splitSequence = BL.CodonSequence(mrnaSequence)
translatedAndAligned = BL.alignseq(splitSequence)
justAminoAcids = BL.translatedSequence(splitSequence)#check this
codonFrequency = BL.codonFreq(splitSequence)# need to edit to incorporate total frequencies
restrictionEnzymeCutSites = BL.restrictionEnzyme('ttgtc', start, end, parsedSequence) #returned as dictionary




#for back end

#totalCodonfrequency = BL.totalCodonFreq(allseq)


