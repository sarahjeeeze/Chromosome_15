import final_version

# put codons in seperate dictionary :) ? & import in to each function or class that uses them 

with open('seq.txt', 'r') as f:
  numbers = f.read().split()

print(numbers)
parsed = final_version.Parse_sequence(numbers)

print(parsed)

gene = final_version.coding_region(53,450, parsed)
print(gene)

intron = final_version.introns(gene, parsed)

print(intron)

introns_stuck = final_version.introns_stuck_together(intron)

print(introns_stuck)


amino_seq = final_version.alignseq(introns_stuck)

print(amino_seq)

codon_freq  = final_version.codon_freq(introns_stuck)

print(codon_freq)


restriction_enzyme = final_version.restriction_enzyme('atcgat',parsed , 59,560)

print(restriction_enzyme)
