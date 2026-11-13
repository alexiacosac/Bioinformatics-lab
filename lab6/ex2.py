
import csv
from pathlib import Path
from urllib.request import urlopen

GENOME_ID = "NC_001604.1"
FASTA_FILE = Path(__file__).parent / "t7_genome.fasta"
CSV_FILE = Path(__file__).parent / "possible_transposons.csv"

MIN_TIR = 6
MAX_TIR = 8
MIN_INTERNAL = 10
MAX_INTERNAL = 150
MIN_TSD = 2
MAX_TSD = 3


def get_genome():
    if FASTA_FILE.exists():
        print("Using existing FASTA file")
        return

    print("Downloading genome from NCBI...")

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={GENOME_ID}&rettype=fasta&retmode=text"
    )

    with urlopen(url, timeout=30) as response:
        content = response.read().decode("utf-8")

    if not content.startswith(">"):
        raise ValueError("Invalid FASTA data")

    FASTA_FILE.write_text(content, encoding="utf-8")

    print("Genome downloaded successfully")


def load_dna():
    sequence = ""

    with open(FASTA_FILE, "r") as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip().upper()

    if not sequence:
        raise ValueError("Empty DNA sequence")

    return sequence


def reverse_complement(dna):
    pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}

    result = ""

    for nucleotide in dna[::-1]:
        result += pairs[nucleotide]

    return result


def search_transposons(dna):
    candidates = {}
    dna_length = len(dna)

    for tir_size in range(MAX_TIR, MIN_TIR - 1, -1):
        print("Checking TIR size:", tir_size)

        last_start = dna_length - 2 * tir_size - MIN_INTERNAL + 1

        for pos in range(last_start):
            left_tir = dna[pos:pos + tir_size]

            if not set(left_tir).issubset(set("ATGC")):
                continue

            right_tir = reverse_complement(left_tir)

            first = pos + tir_size + MIN_INTERNAL
            last = min(
                pos + tir_size + MAX_INTERNAL + tir_size,
                dna_length
            )

            right_pos = dna.find(right_tir, first, last)

            while right_pos != -1:
                element_end = right_pos + tir_size
                internal_size = right_pos - pos - tir_size

                for tsd_size in range(MAX_TSD, MIN_TSD - 1, -1):
                    left_start = pos - tsd_size
                    right_end = element_end + tsd_size

                    if left_start < 0 or right_end > dna_length:
                        continue

                    left_tsd = dna[left_start:pos]
                    right_tsd = dna[element_end:right_end]

                    if left_tsd != right_tsd:
                        continue

                    if not set(left_tsd).issubset(set("ATGC")):
                        continue

                    start = pos + 1
                    end = element_end
                    key = (start, end)

                    candidate = {
                        "start": start,
                        "end": end,
                        "length": end - start + 1,
                        "tir_length": tir_size,
                        "left_tir": left_tir,
                        "right_tir": right_tir,
                        "internal_length": internal_size,
                        "tsd": left_tsd,
                        "tsd_length": tsd_size
                    }

                    previous = candidates.get(key)

                    if previous is None or (
                        tir_size, tsd_size
                    ) > (
                        previous["tir_length"],
                        previous["tsd_length"]
                    ):
                        candidates[key] = candidate

                    break

                right_pos = dna.find(right_tir, right_pos + 1, last)

    return sorted(
        candidates.values(),
        key=lambda item: (
            -item["tir_length"],
            -item["tsd_length"],
            item["start"]
        )
    )


def find_overlaps(candidates):
    overlaps = []
    ordered = sorted(candidates, key=lambda x: x["start"])

    for i, first in enumerate(ordered):
        for second in ordered[i + 1:]:

            if second["start"] > first["end"]:
                break

            if second["end"] <= first["end"]:
                relation = "Nested"
            else:
                relation = "Partial overlap"

            overlaps.append((
                first["start"],
                first["end"],
                second["start"],
                second["end"],
                relation
            ))

    return overlaps


def export_csv(candidates):
    columns = [
        "start", "end", "length",
        "tir_length", "left_tir", "right_tir",
        "internal_length", "tsd", "tsd_length"
    ]

    with open(CSV_FILE, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=columns)

        writer.writeheader()

        for candidate in candidates:
            writer.writerow(candidate)


def main():
    try:
        get_genome()
        dna = load_dna()

    except Exception as error:
        print("Error loading genome:", error)
        print("You can download the FASTA manually from NCBI")
        return

    print("\nGENOME DETAILS")
    print("Accession:", GENOME_ID)
    print("Sequence length:", len(dna), "bp")

    print("\nSEARCHING FOR POSSIBLE TRANSPOSONS")

    candidates = search_transposons(dna)

    print("\nTotal possible transposons:", len(candidates))

    print("\nTOP 20 CANDIDATES")

    for index, candidate in enumerate(candidates[:20], start=1):
        print(
            f"{index}. "
            f"Start: {candidate['start']}, "
            f"End: {candidate['end']}, "
            f"Length: {candidate['length']}, "
            f"TIR: {candidate['left_tir']} / "
            f"{candidate['right_tir']}, "
            f"Internal: {candidate['internal_length']}, "
            f"TSD: {candidate['tsd']}"
        )

    export_csv(candidates)

    print("\nOVERLAPPING CANDIDATES")

    overlaps = find_overlaps(candidates)

    print("Total overlapping pairs:", len(overlaps))

    for overlap in overlaps[:10]:
        print(overlap)

    print("\nResults saved to:", CSV_FILE.name)


if __name__ == "__main__":
    main()
