from pathlib import Path

FASTA_FILE = Path(__file__).parent / "sequence.fasta"


def load_sequence(filename):
    dna = ""

    with open(filename, "r") as file:
        for line in file:
            if not line.startswith(">"):
                dna += line.strip().upper()

    return dna


def check_unit(unit):
    length = len(unit)

    for i in range(1, length):
        if length % i == 0:
            smaller = unit[:i]

            if smaller * (length // i) == unit:
                return False

    return True


def detect_repeats(dna):
    found = []

    for size in [3, 4, 5]:
        position = 0

        while position + size * 2 <= len(dna):
            pattern = dna[position:position + size]

            if not check_unit(pattern):
                position += 1
                continue

            count = 1
            next_position = position + size

            while next_position + size <= len(dna):
                next_pattern = dna[next_position:next_position + size]

                if next_pattern != pattern:
                    break

                count += 1
                next_position += size

            if count >= 2:
                max_count = 10 // size
                used_count = min(count, max_count)

                if used_count >= 2:
                    repeated_dna = pattern * used_count

                    start = position + 1
                    end = position + len(repeated_dna)

                    found.append((
                        start,
                        end,
                        pattern,
                        used_count,
                        repeated_dna
                    ))

                    position = next_position
                    continue

            position += 1

    found.sort()
    return found


def main():
    dna = load_sequence(FASTA_FILE)
    repeats = detect_repeats(dna)

    print("DNA length:", len(dna), "bp")
    print("Number of tandem repeats:", len(repeats))

    for repeat in repeats:
        start, end, pattern, count, sequence = repeat

        print(
            f"{start}-{end}: {sequence} "
            f"(unit={pattern}, copies={count}, "
            f"length={len(sequence)})"
        )


if __name__ == "__main__":
    main()
