import random

with open("lab4/genome.fna", "r") as f:
    dna = "".join(
        row.strip() for row in f
        if not row.startswith(">")
    ).upper()

print("Genome size:", len(dna))

fragments = [
    dna[pos:pos + size]
    for size, pos in (
        (length, random.randrange(len(dna) - length + 1))
        for length in (random.randint(100, 120) for _ in range(2000))
    )
]

print("Number of fragments:", len(fragments))

print("\nRandom DNA fragments:")
for index, fragment in enumerate(fragments[:5], start=1):
    print(f"Fragment {index}: {fragment}")

assembled_dna = "".join(fragments)

print("\nOriginal DNA length:", len(dna))
print("Assembled DNA length:", len(assembled_dna))

print("\nOriginal DNA (first 100 bases):")
print(dna[:100])

print("\nAssembled DNA (first 100 bases):")
print(assembled_dna[:100])

print("\nReconstruction result:")

if dna == assembled_dna:
    print("The original DNA sequence was successfully reconstructed.")
else:
    print("Reconstruction unsuccessful.")
    print("The fragments are unordered and may contain overlaps or missing regions.")