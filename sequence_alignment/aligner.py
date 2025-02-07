from Bio.Align import PairwiseAligner

def align_sequences(seq1, seq2, method="global"):
    aligner = PairwiseAligner()
    aligner.mode = "global" if method == "global" else "local"
    alignments = aligner.align(seq1, seq2)
    return alignments[0]

if __name__ == "__main__":
    seq1 = "GATTACA"
    seq2 = "GCATGCU"
    alignment = align_sequences(seq1, seq2, method="global")
    print(alignment)
