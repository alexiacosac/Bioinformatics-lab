# Exercise 1
# Calculate the melting temperature (Tm) of a DNA sequence.
# Formula: Tm = 4 * (G + C) + 2 * (A + T)
S = input("Enter DNA sequence: ").upper()

A = S.count("A")
T = S.count("T")
G = S.count("G")
C = S.count("C")

Tm = 4 * (G + C) + 2 * (A + T)

print("Melting temperature:", Tm, "°C")