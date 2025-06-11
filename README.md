# Global Sequence Alignment with HMM Emission Profile

This project implements a **global sequence alignment** algorithm using the Needleman-Wunsch technique. It also constructs a simple **profile HMM** by calculating emission probabilities from aligned sequences.

## 🧬 Features

- Global alignment of DNA sequences using match/mismatch/gap penalties.
- Progressive multiple sequence alignment.
- Emission probability computation for each position (HMM Profile).
- Visualization of alignments with match indicators.
- Random DNA sequence generator with realistic biological patterns.
- Alignment score reporting for datasets and random pairs.

## 🧪 Algorithms Used

- **Needleman-Wunsch Algorithm** for global alignment.
- **Progressive alignment** for multiple sequence alignment.
- **Profile HMM** construction by emission frequency estimation.

## 📁 Dataset

- **Dataset A:** 10 sequences used for alignment and display.
- **Dataset B:** 70 sequences used to train the HMM profile.
- **Dataset C:** 20 sequences used to test pairwise alignments.

## 📊 Scoring

- Match: +1
- Mismatch: -1
- Gap: -2

## 📦 Requirements

- Python 3.x
- No external libraries required.

## 🚀 Running the Program

To run the program:

```bash
python3 sequence_alignment.py
