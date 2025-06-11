import random

# Scoring scheme
MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_PENALTY = -2

# Needleman-Wunsch Global Alignment
def global_sequence_alignment(seqA, seqB):
    m, n = len(seqA), len(seqB)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize first row and column with gap penalties
    for i in range(m + 1):
        score_matrix[i][0] = GAP_PENALTY * i
    for j in range(n + 1):
        score_matrix[0][j] = GAP_PENALTY * j

    # Fill in the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_mismatch = MATCH_SCORE if seqA[i - 1] == seqB[j - 1] else MISMATCH_SCORE
            diagonal = score_matrix[i - 1][j - 1] + match_mismatch
            up = score_matrix[i - 1][j] + GAP_PENALTY
            left = score_matrix[i][j - 1] + GAP_PENALTY
            score_matrix[i][j] = max(diagonal, up, left)

    return backtrack(seqA, seqB, score_matrix)

# Backtrack to get aligned sequences
def backtrack(seqA, seqB, score_matrix):
    i, j = len(seqA), len(seqB)
    alignedA, alignedB = [], []

    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + (
                MATCH_SCORE if seqA[i - 1] == seqB[j - 1] else MISMATCH_SCORE):
            alignedA.append(seqA[i - 1])
            alignedB.append(seqB[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + GAP_PENALTY:
            alignedA.append(seqA[i - 1])
            alignedB.append("-")
            i -= 1
        else:
            alignedA.append("-")
            alignedB.append(seqB[j - 1])
            j -= 1

    return "".join(reversed(alignedA)), "".join(reversed(alignedB))

# Align a list of sequences progressively
def multiple_sequence_alignment(sequences):
    aligned_sequences = [sequences[0]]
    for i in range(1, len(sequences)):
        result = global_sequence_alignment(aligned_sequences[i - 1], sequences[i])
        aligned_sequences[i - 1] = result[0]
        aligned_sequences.append(result[1])
    return aligned_sequences

# Pad sequences to the same length
def pad_sequences(aligned_sequences):
    max_length = max(len(seq) for seq in aligned_sequences)
    for i in range(len(aligned_sequences)):
        while len(aligned_sequences[i]) < max_length:
            aligned_sequences[i] += "-"
    return aligned_sequences

# Calculate emission probabilities for HMM profile
def calculate_emission_probabilities(aligned_sequences, position):
    counts = [0] * 5  # A, C, G, T, Gap
    for seq in aligned_sequences:
        symbol = seq[position]
        if symbol == 'A':
            counts[0] += 1
        elif symbol == 'C':
            counts[1] += 1
        elif symbol == 'G':
            counts[2] += 1
        elif symbol == 'T':
            counts[3] += 1
        elif symbol == '-':
            counts[4] += 1
    total = len(aligned_sequences)
    return [count / total for count in counts]

# Initialize and print HMM emission probabilities
def initialize_hmm(aligned_sequences):
    sequence_length = len(aligned_sequences[0])
    for i in range(sequence_length):
        emissions = calculate_emission_probabilities(aligned_sequences, i)
        print(f"Position {i} - Emission Probabilities:")
        print(f"A: {emissions[0]}, C: {emissions[1]}, G: {emissions[2]}, T: {emissions[3]}, Gap: {emissions[4]}")

# Train HMM model from a set of aligned sequences
def train_hmm(dataset_b):
    aligned_sequences = multiple_sequence_alignment(dataset_b)
    aligned_sequences = pad_sequences(aligned_sequences)
    initialize_hmm(aligned_sequences)

# Score alignment and print with visual indicator
def print_alignment_scores_and_paths(dataset):
    for i in range(len(dataset)):
        for j in range(i + 1, len(dataset)):
            seqA, seqB = dataset[i], dataset[j]
            result = global_sequence_alignment(seqA, seqB)
            score = score_alignment(seqA, seqB, result[0], result[1])

            print(f"Alignment score between Seq{i} and Seq{j}: {score}")
            print("Aligned Sequences:")
            line = ''.join(['|' if a == b and a != '-' else ' ' for a, b in zip(result[0], result[1])])
            print(result[0])
            print(line)
            print(result[1])

# Compute total alignment score
def score_alignment(seqA, seqB, alignedA, alignedB):
    score = 0
    for a, b in zip(alignedA, alignedB):
        if a == b:
            score += MATCH_SCORE
        elif a == '-' or b == '-':
            score += GAP_PENALTY
        else:
            score += MISMATCH_SCORE
    return score

# Generate a random sequence of A, C, G, T
def generate_random_sequence():
    symbols = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(symbols) for _ in range(random.randint(5, 30)))

# Main method to drive the pipeline
def main():
    all_sequences = []
    symbols = ['A', 'C', 'G', 'T']
    patterns = ["ATTAGA", "ACGCATTT", "AGGACTCAA", "ATTTCAGT"]

    # Generate dataset with mutations
    for _ in range(100):
        sequence_parts = []

        prefix = ''.join(random.choice(symbols) for _ in range(random.randint(1, 3)))
        sequence_parts.append(prefix)

        for pattern in patterns:
            symbol_set = list(set(pattern))
            mod_type = random.randint(0, 2)

            if mod_type == 0:
                sequence_parts.append(pattern)
            elif mod_type == 1:
                idx = random.randint(0, len(pattern) - 1)
                if random.choice([0, 1]) == 0:
                    pattern = pattern[:idx] + pattern[idx + 1:]  # deletion
                else:
                    new_char = random.choice(symbol_set)
                    while new_char == pattern[idx]:
                        new_char = random.choice(symbol_set)
                    pattern = pattern[:idx] + new_char + pattern[idx + 1:]
                sequence_parts.append(pattern)
            else:
                idx1, idx2 = sorted(random.sample(range(len(pattern)), 2))
                deleted = 0
                for idx in [idx1, idx2]:
                    if random.choice([0, 1]) == 0 and len(pattern) > 1:
                        pattern = pattern[:idx - deleted] + pattern[idx - deleted + 1:]
                        deleted += 1
                    else:
                        new_char = random.choice(symbol_set)
                        while new_char == pattern[idx - deleted]:
                            new_char = random.choice(symbol_set)
                        pattern = pattern[:idx - deleted] + new_char + pattern[idx - deleted + 1:]
                sequence_parts.append(pattern)

        suffix = ''.join(random.choice(symbols) for _ in range(random.randint(1, 2)))
        sequence_parts.append(suffix)

        all_sequences.append(''.join(sequence_parts))

    random.shuffle(all_sequences)
    dataset_a = all_sequences[:10]
    dataset_b = all_sequences[10:80]
    dataset_c = all_sequences[80:]

    print("\n--- Dataset A ---")
    for s in dataset_a:
        print(s)

    print("\n--- Dataset B ---")
    for s in dataset_b:
        print(s)

    print("\n--- Dataset C ---")
    for s in dataset_c:
        print(s)

    print("\n--- Aligned Sequences in Dataset A ---")
    aligned = pad_sequences(multiple_sequence_alignment(dataset_a))
    for seq in aligned:
        print(seq)

    print("\n--- Emission Probabilities using HMM Profile ---")
    train_hmm(dataset_b)

    print("\n--- Alignment Scores and Paths for Dataset C ---")
    print_alignment_scores_and_paths(dataset_c)

    print("\n--- Alignment Scores for 20 Random Sequence Pairs ---")
    for _ in range(20):
        seq_a = generate_random_sequence()
        seq_b = generate_random_sequence()
        result = global_sequence_alignment(seq_a, seq_b)
        score = score_alignment(seq_a, seq_b, result[0], result[1])
        print(f"Alignment score: {score}")
        print("Aligned Sequences:")
        print(result[0])
        print(result[1])

if __name__ == "__main__":
    main()
