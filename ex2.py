# Exercise 2: Calculate the melting temperature of a DNA sequence
# Formula: Tm = 81.5 + 16.6 * log10(Na) + 41 * GC - 600 / length
import math

S = input("Enter DNA sequence: ").upper()
Na = float(input("Enter Na+ concentration: "))

length = len(S)
GC = (S.count("G") + S.count("C")) / length

Tm = 81.5 + 16.6 * math.log10(Na) + 41 * GC - 600 / length

print("Melting temperature:", round(Tm, 2), "°C")