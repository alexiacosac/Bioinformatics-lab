import argparse
from pathlib import Path
from statistics import mean

import matplotlib.pyplot as plt


WINDOW_SIZE = 30
TEST_SEQUENCE = (
    "CGGACTGATCTATCT"
    + "A" * 27
    + "CGTAGCATCTATCGATCTATCTAGCGATCTA"
)


def read_fasta(path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    headers = [line for line in lines if line.startswith(">")]

    if len(headers) != 1:
        raise ValueError("The FASTA file must contain one promoter sequence.")

    sequence = "".join(
        line.strip() for line in lines if not line.startswith(">")
    ).upper()

    return headers[0][1:].split()[0], sequence


def gc_percent(sequence):
    return 100 * (sequence.count("C") + sequence.count("G")) / len(sequence)


def kappa_ic(window):
    n = len(window)
    coincidence_rates = []

    for shift in range(1, n):
        matches = sum(
            window[i] == window[i + shift]
            for i in range(n - shift)
        )
        coincidence_rates.append(matches / (n - shift))

    return 100 * mean(coincidence_rates)


def calculate_pattern(sequence, relative_gc):
    if len(sequence) < WINDOW_SIZE:
        raise ValueError(f"The sequence must have at least {WINDOW_SIZE} bases.")

    if set(sequence) - set("ACGT"):
        raise ValueError("The sequence must contain only A, C, G and T.")

    total_gc = gc_percent(sequence)
    x_values = []
    y_values = []

    for start in range(len(sequence) - WINDOW_SIZE + 1):
        window = sequence[start:start + WINDOW_SIZE]
        x = gc_percent(window)

        if relative_gc:
            x = x * total_gc / 100

        x_values.append(x)
        y_values.append(kappa_ic(window))

    return x_values, y_values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, help="FASTA file with one promoter")
    parser.add_argument("--relative-gc", action="store_true")
    args = parser.parse_args()

    sequences = [("Artificial test", TEST_SEQUENCE)]
    if args.fasta:
        sequences.append(read_fasta(args.fasta))

    x_label = "Relative C+G" if args.relative_gc else "C+G (%)"
    pattern_figure, pattern_axis = plt.subplots(layout="constrained")
    center_figure, center_axis = plt.subplots(layout="constrained")

    for name, sequence in sequences:
        x_values, y_values = calculate_pattern(sequence, args.relative_gc)
        center_x = mean(x_values)
        center_y = mean(y_values)

        print(f"\n{name}")
        print(f"Sequence length: {len(sequence)}")
        print(f"Windows: {len(x_values)}")
        print(f"Whole sequence GC: {gc_percent(sequence):.2f}%")
        print(f"First window: X = {x_values[0]:.2f}, IC = {y_values[0]:.2f}")
        print(f"Center: ({center_x:.2f}, {center_y:.2f})")

        pattern_axis.scatter(x_values, y_values, s=18, alpha=0.65, label=name)
        center_axis.scatter(center_x, center_y, s=95, marker="X", label=name)

    for axis, title in [
        (pattern_axis, "DNA patterns"),
        (center_axis, "Pattern centers"),
    ]:
        axis.set(title=title, xlabel=x_label, ylabel="Kappa IC (%)")
        axis.grid(alpha=0.25)
        axis.legend()

    pattern_figure.savefig("dna_patterns.png", dpi=180)
    center_figure.savefig("pattern_centers.png", dpi=180)


if __name__ == "__main__":
    main()