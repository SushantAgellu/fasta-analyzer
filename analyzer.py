from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt

def analyze_fasta(file_path):
    total_sequences = 0
    total_length = 0
    gc_values = []

    for record in SeqIO.parse(file_path, "fasta"):
        total_sequences += 1
        seq_len = len(record.seq)
        total_length += seq_len
        gc_content = ((record.seq.count("G") + record.seq.count("C")) / seq_len) * 100
        gc_values.append((record.id, gc_content))

        # --- Codon usage ---
        codons = [str(record.seq[i:i+3]) for i in range(0, len(record.seq) - 2, 3)]
        codon_count = Counter(codons)
        top_codons = codon_count.most_common(5)

        # --- ORF detection ---
        start_codon = "ATG"
        stops = {"TAA", "TAG", "TGA"}
        orfs = []
        for i in range(len(record.seq) - 2):
            codon = record.seq[i:i+3]
            if codon == start_codon:
                for j in range(i + 3, len(record.seq) - 2, 3):
                    stop = record.seq[j:j+3]
                    if stop in stops:
                        orfs.append((i + 1, j + 3, j - i))
                        break

        print(f"\n{record.id}:")
        print(f"  Length: {seq_len}")
        print(f"  GC Content: {gc_content:.2f}%")
        print(f"  Top 5 Codons: {top_codons}")
        print(f"  ORFs Found: {len(orfs)}")

    print("\nSummary:")
    print(f"  Total sequences: {total_sequences}")
    print(f"  Total length: {total_length}")

    # --- Plot GC content ---
    plot_gc(gc_values)

def plot_gc(gc_values):
    ids = [x[0] for x in gc_values]
    gcs = [x[1] for x in gc_values]
    plt.bar(ids, gcs)
    plt.xlabel("Sequence ID")
    plt.ylabel("GC%")
    plt.title("GC Content per Sequence")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    path = input("Enter FASTA file path: ").strip()
    analyze_fasta(path)
