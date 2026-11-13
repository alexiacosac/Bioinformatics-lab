
import random
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path

FASTA_FILE = Path(__file__).parent / "sequence.fasta"


def load_dna(file_name):
    sequence = ""

    with open(file_name, "r") as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip().upper()

    return sequence


def generate_fragments(sequence, number=10):
    fragments = []

    for i in range(number):
        size = random.randint(100, min(3000, len(sequence)))
        position = random.randint(0, len(sequence) - size)

        fragment = sequence[position:position + size]
        fragments.append(fragment)

    return fragments


def calculate_distance(size):
    max_size = 3000
    min_size = 100

    distance = 1 + 8 * (
        math.log10(max_size / size) /
        math.log10(max_size / min_size)
    )

    return distance


def simulate_gel(sizes, distances):
    fig, ax = plt.subplots(figsize=(12, 8))

    ax.set_facecolor("black")
    ax.set_xlim(-0.6, 10.6)
    ax.set_ylim(10, -0.5)

    ax.set_title("DNA Gel Electrophoresis")
    ax.set_ylabel("Migration distance (cm)")
    ax.set_xticks(range(11))
    ax.set_xticklabels(["L"] + list(range(1, 11)))

    ladder_sizes = [3000, 2500, 2000, 1500, 1000,
                    750, 500, 250, 100]

    for size in ladder_sizes:
        distance = calculate_distance(size)

        ax.plot(
            [-0.35, 0.35],
            [distance, distance],
            color="cyan",
            linewidth=5
        )

        ax.text(
            -0.42, distance, f"{size}",
            color="white",
            ha="right",
            va="center",
            fontsize=9
        )

    bands = []
    labels = []

    for i in range(10):
        band, = ax.plot(
            [i + 0.65, i + 1.35],
            [0, 0],
            color="orange",
            linewidth=5
        )
        bands.append(band)

        label = ax.text(
            i + 1, 0.25, str(sizes[i]),
            color="white",
            ha="center",
            fontsize=9
        )
        labels.append(label)

    def animate(frame):
        fraction = frame / 50

        for i in range(10):
            current_y = distances[i] * fraction

            bands[i].set_ydata([current_y, current_y])
            labels[i].set_position((i + 1, current_y + 0.25))

        return bands + labels

    animation = FuncAnimation(
        fig,
        animate,
        frames=51,
        interval=50,
        repeat=False
    )

    plt.show()

    animate(50)
    fig.savefig("gel_electrophoresis.png", dpi=150,
                bbox_inches="tight")


def main():
    dna = load_dna(FASTA_FILE)

    fragments = generate_fragments(dna)

    sizes = []
    distances = []

    print("DNA sequence length:", len(dna), "bp")

    for index, fragment in enumerate(fragments):
        size = len(fragment)
        distance = calculate_distance(size)

        sizes.append(size)
        distances.append(distance)

        print(
            f"Fragment {index + 1}: {size} bp | "
            f"Distance: {distance:.2f} cm"
        )

    simulate_gel(sizes, distances)


if __name__ == "__main__":
    main()
