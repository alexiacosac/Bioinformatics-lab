
import random

DNA_LENGTH = 301


def generate_dna(length):
    bases = "ACGT"
    dna = ""

    for i in range(length):
        dna += random.choice(bases)

    return dna


def get_reverse_complement(sequence):
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}

    result = ""
    for base in sequence[::-1]:
        result += pairs[base]

    return result


def insert_transposons(length, positions, types):
    dna = list(generate_dna(length))
    used_positions = set()

    for name, tir, tsd in types:
        start, end = positions[name]

        first_marker = tsd + tir
        last_marker = get_reverse_complement(tir) + tsd

        if end - start + 1 < len(first_marker) + len(last_marker):
            raise ValueError("Transposon length is too small")

        left_index = start - 1
        right_index = end - len(last_marker)

        for index, marker in [(left_index, first_marker),
                              (right_index, last_marker)]:

            indexes = set(range(index, index + len(marker)))

            if used_positions.intersection(indexes):
                raise ValueError("Transposon markers overlap")

            used_positions.update(indexes)
            dna[index:index + len(marker)] = list(marker)

    return "".join(dna)


def find_transposons(dna, types):
    detected = []

    for name, tir, tsd in types:
        tir_end = get_reverse_complement(tir)
        tir_length = len(tir)
        tsd_length = len(tsd)

        for i in range(len(dna)):
            if dna[i:i + tir_length] != tir:
                continue

            start = i - tsd_length

            if start < 0:
                continue

            if dna[start:i] != tsd:
                continue

            for j in range(i + tir_length, len(dna)):
                if dna[j:j + tir_length] != tir_end:
                    continue

                end = j + tir_length + tsd_length

                if dna[j + tir_length:end] == tsd:
                    detected.append((name, start + 1, end))

    detected = list(set(detected))
    detected.sort(key=lambda item: item[1])

    return detected


def analyze_overlaps(transposons):
    overlaps = 0

    for i in range(len(transposons)):
        name1, start1, end1 = transposons[i]

        for j in range(i + 1, len(transposons)):
            name2, start2, end2 = transposons[j]

            if start1 <= end2 and start2 <= end1:
                overlaps += 1

                if start1 <= start2 and end1 >= end2:
                    print(f"{name2} is inside {name1}")

                elif start2 <= start1 and end2 >= end1:
                    print(f"{name1} is inside {name2}")

                else:
                    print(f"{name1} overlaps with {name2}")

    if overlaps == 0:
        print("No overlapping transposons")


def run_case(title, positions, types):
    print("\n" + "=" * 60)
    print("CASE:", title)
    print("=" * 60)

    dna = insert_transposons(DNA_LENGTH, positions, types)

    print("\nGENERATED DNA:")
    print(dna)
    print("DNA length:", len(dna))

    detected = find_transposons(dna, types)

    print("\nDETECTED TRANSPOSONS:")

    for name, start, end in detected:
        length = end - start + 1

        print(f"{name}: start={start}, end={end}, length={length}")

    print("\nOVERLAP ANALYSIS:")
    analyze_overlaps(detected)

    print("\nTRANSPOSON SEQUENCES:")

    for name, start, end in detected:
        fragment = dna[start - 1:end]

        print(f"{name} [{start}-{end}]")
        print(fragment)
        print()

    expected = positions.copy()
    actual = {name: (start, end) for name, start, end in detected}

    if actual == expected:
        print("RESULT: All transposons detected correctly")
    else:
        print("RESULT: Detected positions do not match")


def main():
    random.seed(42)

    types = [
        ("T1", "ATGCCGA", "TAC"),
        ("T2", "CGTATGC", "GGA"),
        ("T3", "GCTACGA", "ACT"),
        ("T4", "TGCAGAT", "CCA")
    ]

    scenarios = {
        "1. SEPARATED": {
            "T1": (15, 55),
            "T2": (80, 120),
            "T3": (150, 190),
            "T4": (225, 265)
        },

        "2. OVERLAPPING": {
            "T1": (30, 110),
            "T2": (85, 165),
            "T3": (185, 225),
            "T4": (250, 290)
        },

        "3. NESTED": {
            "T1": (30, 185),
            "T2": (65, 155),
            "T3": (95, 130),
            "T4": (225, 265)
        },

        "4. THREE INTERSECTING": {
            "T1": (25, 120),
            "T2": (60, 155),
            "T3": (95, 190),
            "T4": (230, 270)
        }
    }

    for title, positions in scenarios.items():
        run_case(title, positions, types)


if __name__ == "__main__":
    main()
