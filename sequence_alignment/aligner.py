from Bio.Align import PairwiseAligner

def align_sequences(seq1, seq2, method="global"):
    aligner = PairwiseAligner()
    aligner.mode = "global" if method == "global" else "local"
    return aligner.align(seq1, seq2)[0]

if __name__ == "__main__":
    s1 = "ACGTAG"
    s2 = "ACGTG"
    print(align_sequences(s1, s2, method="global"))
