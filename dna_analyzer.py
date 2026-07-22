CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}  #CODON_TABLE: maps every possible 3-letter DNA codon to the amino acid it codes for.


def parse_fasta(filepath):
    records = []
    header = None
    seq_lines = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, ''.join(seq_lines)))
    return records
    # Reads a FASTA file and splits it into records.
# Each record is a (header, sequence) pair - the header is the line after '>',
# and the sequence is all DNA letters until the next '>' (joined into one string).



def count_bases(seq):
    counts = {}
    for base in seq:
        if base in counts:
            counts[base] = counts[base] + 1
        else:
            counts[base] = 1
    return counts
    # Counts how many times each base (A, T, G, C) appears in the sequence.



def gc_content(seq):
    counts = count_bases(seq)
    g = counts.get('G', 0)
    c = counts.get('C', 0)
    total = len(seq)
    percentage = (g + c) / total * 100
    return percentage
    # Calculates what percentage of the sequence is made up of G and C bases.



def reverse_complement(seq):
    partner = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    complemented = ''
    for base in seq:
        complemented = complemented + partner[base]
    reversed_seq = complemented[::-1]
    return reversed_seq
#Builds the "opposite strand" of the DNA: swaps each base for its pair
# (A<->T, G<->C) and then reverses the whole sequence.



def transcribe(seq):
    mrna = seq.replace('T', 'U')
    return mrna
# Converts DNA into mRNA by replacing every 'T' with 'U'.



def translate(seq):
    protein = ''
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        amino_acid = CODON_TABLE[codon]
        if amino_acid == '*':
            break
        protein = protein + amino_acid
    return protein
# Converts the sequence into a protein by reading it 3 letters (1 codon) at a time
# and looking each codon up in CODON_TABLE. Stops as soon as a stop codon ('*') is found.



def find_motif(seq, motif):
    positions = []
    start = 0
    while True:
        index = seq.find(motif, start)
        if index == -1:
            break
        positions.append(index)
        start = index + 1
    return positions


if __name__ == "__main__":
    FILE_PATH = "sequences.fasta"   # path to your FASTA file
    MOTIF = "ATG"                   # set to None if you don't want to search for a motif
    SHOW_TRANSLATION = True         # set to False if you don't want the protein translation

    records = parse_fasta(FILE_PATH)
    for header, seq in records:
        print("\n===", header, "===")
        print("Length:", len(seq), "bp")
        print("Base counts:", count_bases(seq))
        print("GC content:", gc_content(seq), "%")
        print("Reverse complement:", reverse_complement(seq))
        print("mRNA:", transcribe(seq))
        if SHOW_TRANSLATION:
            print("Protein:", translate(seq))
        if MOTIF:
            positions = find_motif(seq, MOTIF)
            print(f"Motif '{MOTIF}' found at:", positions)
