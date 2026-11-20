from pathlib import Path
from urllib.request import urlopen

import matplotlib.pyplot as plt


GENOME_ID = "NC_001416.1"
FASTA_FILE = Path(__file__).resolve().parent / "sequence.fasta"

# Each tuple contains the recognition sequence and the number
# of bases before the cut on the 5' -> 3' strand.
ENZYMES = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
}


def read_genome():
    if FASTA_FILE.exists():
        fasta = FASTA_FILE.read_text(encoding="utf-8")
    else:
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=nuccore&id={GENOME_ID}&rettype=fasta&retmode=text"
        )
        print(f"Downloading {GENOME_ID} from NCBI...")
        with urlopen(url, timeout=30) as response:
            fasta = response.read().decode("utf-8")

    lines = fasta.splitlines()
    if not lines or not lines[0].startswith(f">{GENOME_ID} "):
        raise ValueError(f"Expected a FASTA record for {GENOME_ID}.")

    dna = "".join(line.strip() for line in lines[1:]).upper()
    if len(dna) != 48502 or set(dna) - set("ACGT"):
        raise ValueError("The lambda genome sequence is incomplete or invalid.")

    if not FASTA_FILE.exists():
        FASTA_FILE.write_text(fasta, encoding="utf-8")

    return dna


def find_cuts(dna, recognition_sequence, cut_offset):
    cuts = []
    start = dna.find(recognition_sequence)

    while start != -1:
        cuts.append(start + cut_offset)
        start = dna.find(recognition_sequence, start + 1)

    return cuts


def make_fragments(dna_length, cuts):
    boundaries = [0, *cuts, dna_length]
    return [
        (left + 1, right, right - left)
        for left, right in zip(boundaries, boundaries[1:])
    ]


def print_digest(name, cuts, fragments):
    print(f"\n{name}")
    print(f"Number of cleavages: {len(cuts)}")
    print("Cuts after bases:", ", ".join(map(str, cuts)) or "none")
    print(f"Number of fragments: {len(fragments)}")

    for index, (start, end, length) in enumerate(fragments, start=1):
        print(f"  Fragment {index:2}: {start:5}–{end:5} | {length:5} bp")


def plot_gel(lanes):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_facecolor("#14213d")

    for lane, fragments in enumerate(lanes.values(), start=1):
        for _, _, length in fragments:
            ax.hlines(length, lane - 0.28, lane + 0.28,
                      color="#79e7ff", linewidth=3)

    ax.set_yscale("log")
    ax.set_ylim(200, 65000)
    ax.set_xlim(0.5, len(lanes) + 0.5)
    ax.set_xticks(range(1, len(lanes) + 1), list(lanes))
    ax.xaxis.tick_top()
    ax.set_ylabel("Fragment size (bp, log scale)")
    ax.set_title("Simulated agarose gel — bacteriophage lambda", pad=25)

    fig.tight_layout()
    fig.savefig("restriction_gel.png", dpi=200)
    plt.show()


def main():
    dna = read_genome()
    print(f"Genome: {GENOME_ID} | Length: {len(dna)} bp")

    cuts_by_enzyme = {
        name: find_cuts(dna, sequence, offset)
        for name, (sequence, offset) in ENZYMES.items()
    }
    combined_cuts = sorted({
        cut for cuts in cuts_by_enzyme.values() for cut in cuts
    })

    digests = {
        "Uncut": [],
        **cuts_by_enzyme,
        "EcoRI + BamHI + HindIII": combined_cuts,
    }

    lanes = {}
    for name, cuts in digests.items():
        fragments = make_fragments(len(dna), cuts)
        print_digest(name, cuts, fragments)
        lanes[name] = fragments

    plot_gel(lanes)


if __name__ == "__main__":
    main()