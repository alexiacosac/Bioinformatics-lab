import random


DNA_LENGTH = 320
REPEAT_LENGTH = 12
MIN_GAP = 10

# Inclusive coordinates, starting at 1.
# T1 and T3 overlap. T2 is nested inside T1.
# T4 and T5 are separate and occur one after the other.
TRANSPOSONS = [
    ("T1", 20, 125),
    ("T2", 40, 75),
    ("T3", 85, 180),
    ("T4", 195, 240),
    ("T5", 255, 300),
]

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def generate_dna():
    rng = random.Random(42)
    dna = rng.choices("ACGT", k=DNA_LENGTH)

    for _, start, end in TRANSPOSONS:
        left_start = start - 1
        right_start = end - REPEAT_LENGTH

        left_repeat = "".join(rng.choices("ACGT", k=REPEAT_LENGTH))
        right_repeat = reverse_complement(left_repeat)

        dna[left_start:left_start + REPEAT_LENGTH] = left_repeat
        dna[right_start:right_start + REPEAT_LENGTH] = right_repeat

    return "".join(dna)


def detect_transposons(dna):
    detected = []
    last_repeat_start = len(dna) - REPEAT_LENGTH

    for left_start in range(last_repeat_start + 1):
        left_end = left_start + REPEAT_LENGTH
        left_repeat = dna[left_start:left_end]
        expected_right_repeat = reverse_complement(left_repeat)

        # Leave at least MIN_GAP bases between the two repeats.
        first_right_start = left_end + MIN_GAP

        for right_start in range(first_right_start, last_repeat_start + 1):
            right_end = right_start + REPEAT_LENGTH
            right_repeat = dna[right_start:right_end]

            if right_repeat == expected_right_repeat:
                start = left_start + 1
                end = right_end
                length = end - start + 1
                detected.append((start, end, length))

    return detected


def main():
    dna = generate_dna()
    print(f"Generated DNA ({len(dna)} bp):")
    for start in range(0, len(dna), 80):
        print(dna[start:start + 80])

    print("\nDetected transposons:")
    print(f"{'Start':>6} {'End':>6} {'Length (bp)':>12}")

    for start, end, length in detect_transposons(dna):
        print(f"{start:>6} {end:>6} {length:>12}")


if __name__ == "__main__":
    main()