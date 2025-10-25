from Bio import SeqIO

def analyze_fasta(file_path):
    total_sequences = 0
    total_length = 0

    for record in SeqIO.parse(file_path, "fasta"):
        total_sequences += 1
        seq_len = len(record.seq)
        total_length += seq_len
        gc_content = ((record.seq.count("G") + record.seq.count("C")) / seq_len )* 100
        print(f"{record.id}: length={seq_len}, GC={gc_content:.2f}%")

    print("\nSummary:")
    print(f"Total sequences: {total_sequences}")
    print(f"Total length: {total_length}")

if __name__ == "__main__":
    path = input("Enter FASTA file path: ").strip()
    analyze_fasta(path)
