def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    scores = [[0] * cols for _ in range(rows)]

    for i in range(rows):
        scores[i][0] = i * gap

    for j in range(cols):
        scores[0][j] = j * gap

    for i in range(1, rows):
        for j in range(1, cols):
            pair_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal = scores[i - 1][j - 1] + pair_score
            up = scores[i - 1][j] + gap
            left = scores[i][j - 1] + gap
            scores[i][j] = max(diagonal, up, left)

    aligned1 = []
    aligned2 = []
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            pair_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            if scores[i][j] == scores[i - 1][j - 1] + pair_score:
                aligned1.append(seq1[i - 1])
                aligned2.append(seq2[j - 1])
                i -= 1
                j -= 1
                continue

        if i > 0 and scores[i][j] == scores[i - 1][j] + gap:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    return "".join(reversed(aligned1)), "".join(reversed(aligned2)), scores


if __name__ == "__main__":
    seq1 = "ACCGTGAAGCCAATAC"
    seq2 = "AGCGTGCAGCCAATAC"

    aligned1, aligned2, scores = needleman_wunsch(seq1, seq2)

    print("Score matrix:")
    print("        " + "".join(f"{base:>4}" for base in seq2))
    for label, row in zip("-" + seq1, scores):
        print(f"{label:>4}" + "".join(f"{score:>4}" for score in row))

    markers = "".join(
        " " if "-" in (a, b) else "|" if a == b else "."
        for a, b in zip(aligned1, aligned2)
    )
    print("\nOptimal alignment:")
    print(aligned1)
    print(markers)
    print(aligned2)
    print(f"Alignment score: {scores[-1][-1]}")