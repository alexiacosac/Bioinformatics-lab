from collections import Counter
import matplotlib.pyplot as plt

filename = "GCF_001482955.1_ViralProj306892_genomic.fna"

with open(filename, "r") as file:
    genome = ""
    for line in file:
        if not line.startswith(">"):
            genome += line.strip().upper()

codons = []

for i in range(len(genome) - 2):
    codon = genome[i:i+3]

    if all(letter in "ATGC" for letter in codon):
        codons.append(codon)

frequency = Counter(codons)

top10 = frequency.most_common(10)

print("Top 10 most frequent codons:")

for codon, count in top10:
    print(codon, ":", count)

names = [codon for codon, count in top10]
counts = [count for codon, count in top10]

plt.bar(names, counts)
plt.xlabel("Codons")
plt.ylabel("Frequency")
plt.title("Top 10 Most Frequent Codons")
plt.show()
