#!/bin/python3

import argparse
import pysam
from Bio.Seq import Seq

# Funktion, um das Reverse-Complement zu erzeugen
def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())

# Lesen der FASTA-Datei und Speichern der Header und Sequenzen
def read_fasta(fasta_file):
    fasta_data = []
    header = None
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]  # Header ohne das ">"-Zeichen
            else:
                sequence = line.strip()
                rev_comp = reverse_complement(sequence)
                fasta_data.append((header, sequence, rev_comp))
    return fasta_data

# Zählen der Vorkommen einer Sequenz in der BAM-Datei (Teilstring-Vergleich)
def count_sequences_in_bam(bam_file, fasta_data):
    sequence_counts = {header: {'forward_count': 0, 'reverse_count': 0} for header, _, _ in fasta_data}
    
    # Öffnen und Parsen der BAM-Datei
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            read_sequence = read.query_sequence
            # Vergleichen mit jeder gespeicherten FASTA-Sequenz (Teilstring-Vergleich)
            for header, sequence, rev_comp in fasta_data:
                if sequence in read_sequence:
                    sequence_counts[header]['forward_count'] += 1
                if rev_comp in read_sequence:
                    sequence_counts[header]['reverse_count'] += 1
    
    return sequence_counts

# Ausgabe der Ergebnisse
def output_results(sequence_counts, sample_id):
    print("ID;Hybrid;fwd_count;rev_count;sum_count")
    for header, counts in sequence_counts.items():
        forward_count = counts['forward_count']
        reverse_count = counts['reverse_count']
        total_count = forward_count + reverse_count
        print(f"{sample_id};{header};{forward_count};{reverse_count};{total_count}")

# Hauptfunktion
def main(fasta_file, bam_file, sample_id):
    fasta_data = read_fasta(fasta_file)
    sequence_counts = count_sequences_in_bam(bam_file, fasta_data)
    output_results(sequence_counts, sample_id)

if __name__ == "__main__":
    # Kommandozeilenargumente parsen
    parser = argparse.ArgumentParser(description="BAM and FASTA Sequence Parser")
    parser.add_argument('-f', '--fasta', required=True, help="Pfad zur FASTA-Datei")
    parser.add_argument('-b', '--bam', required=True, help="Pfad zur BAM-Datei")
    parser.add_argument('-i', '--id', required=True, help="Sample ID für die Ausgabe")

    args = parser.parse_args()

    # Aufrufen der Hauptfunktion mit den übergebenen Argumenten
    main(args.fasta, args.bam, args.id)

