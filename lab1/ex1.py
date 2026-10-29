seq = input("Sequence:  ")

for character in set(seq):
    freq = seq.count(character) / 20
    print(character, ":", freq)