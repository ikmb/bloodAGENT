#!/bin/python3

import argparse
import pyBigWig
import pandas as pd
import sys
import os

# Funktion, um Coverage aus BigWig zu berechnen
def calculate_coverage(bigwig_path, regions):
    bw = pyBigWig.open(bigwig_path)
    coverages = []
    for _, row in regions.iterrows():
        chrom = row['Chromosome']
        start = row['Start']
        end = row['End']
        # Durchschnittliche Coverage des Bereichs berechnen
        coverage = bw.stats(chrom, start, end, type="mean")[0]
        if coverage is not None:
            coverages.append(coverage)
        else:
            coverages.append(0)
    bw.close()
    return coverages

def main():
    # Argumente einlesen
    parser = argparse.ArgumentParser(description="Calculate coverage and CN state for RHCE exon 2.")
    parser.add_argument("-c", "--coverage", required=True, help="Path to BigWig file.")
    parser.add_argument("-b", "--build", required=True, help="genome build (hg38 or hg19)")
    parser.add_argument("-p", "--pipeline", required=True, help="What kind of secondary analysis pipeline is used? (GATK, TGS, Dragen")
    parser.add_argument("-o", "--out", required=True, help="Path to output VCF file.")
    parser.add_argument("-a", "--annotation", required=True, help="Path to annotation file.")

    args = parser.parse_args()

    if args.build not in ("hg19", "hg38"):
        sys.exit("Error: Unsupported build. Please specify 'hg19' or 'hg38'.")
    if not os.path.exists(args.coverage):
        sys.exit(f"Error: The file '{args.coverage}' does not exist.")
    if not os.path.exists(args.annotation):
        sys.exit(f"Error: The file '{args.annotation}' does not exist.")
    try:
        with open(args.out, 'w') as f:
            pass  # Testweise die Datei öffnen und sofort schließen
    except IOError:
        sys.exit(f"Error: The output file '{args.out}' could not be created.")
        
    REF_Allele = "G"
    ALT_Allele = "GAGCTATGATTGTACCACTGGGAAGTGACAAAGGGCACCCTGGGGGATTTCAAATGGTGGTGGCCCTGGTTTGGTGTTGCTGCCAGGTGAGTCCTTAAGCTATAGCAA"
    POS_Allele = 25405596

    match args.pipeline:
        case "GATK":
            REF_Allele = "G"
            ALT_Allele = "GAGCTATGATTGTACCACTGGGAAGTGACAAAGGGCACCCTGGGGGATTTCAAATGGTGGTGGCCCTGGTTTGGTGTTGCTGCCAGGTGAGTCCTTAAGCTATAGCAA"
            POS_Allele = 25405596
            if args.build == "hg19":
                POS_Allele = 25732079
        case "TGS":
            REF_Allele = "T"
            ALT_Allele = "TGCAATGAGCTATGATTGTACCACTGGGAAGTGACAAAGGGCACCCTGGGGGATTTCAAATGGTGGTGGCCCTGGTTTGGTGTTGCTGCCAGGTGAGTCCTTAAGCTATA"
            POS_Allele = 25405592
            if args.build == "hg19":
                POS_Allele = 25732083
        case "Dragen":
            REF_Allele = "G"
            ALT_Allele = "GAGCTATGATTGTACCACTGGGAAGTGACAAAGGGCACCCTGGGGGATTTCAAATGGTGGTGGCCCTGGTTTGGTGTTGCTGCCAGGTGAGTCCTTAAGCTATAGCAA"
            POS_Allele = 25405596
            if args.build == "hg19":
                POS_Allele = 25732079
        case _:
            sys.exit("Unknown pipeline. Please specify 'GATK', 'TGS', or 'Dragen'.")
          
    # Annotation-Datei laden
    annotation = pd.read_csv(args.annotation, sep="\t")

    # RHCE-Annotation extrahieren
    rhce_exons = annotation[annotation['System'] == 'RHCE']
    rhce_exons = rhce_exons[['refGene.chrom', 'refGene.exonStarts', 'refGene.exonEnds']]
    
    # Exon-Koordinaten verarbeiten
    regions = []
    for _, row in rhce_exons.iterrows():
        chrom = row['refGene.chrom']
        starts = list(map(int, row['refGene.exonStarts'].strip(",").split(",")))
        ends = list(map(int, row['refGene.exonEnds'].strip(",").split(",")))
        for start, end in zip(starts, ends):
            regions.append({'Chromosome': chrom, 'Start': start, 'End': end})

    regions_df = pd.DataFrame(regions)

    # Coverage der RHCE-Exons berechnen
    regions_df['Coverage'] = calculate_coverage(args.coverage, regions_df)

    # Durchschnittliche Coverage aller RHCE-Exons
    mean_rhce_coverage = regions_df['Coverage'].mean()

    # Coverage von Exon 2 (Exon 2 ist der achte Eintrag in der Liste)
    exon_2 = regions_df.iloc[8]
    mean_exon_2_coverage = exon_2['Coverage']

    # CNstate berechnen
    cn_state = (mean_exon_2_coverage / mean_rhce_coverage) * 2

    #print(f"Mean RHCE Coverage: {mean_rhce_coverage:.2f}")
    #print(regions_df.to_csv(index=False, sep='\t'))
    # Ergebnis in VCF speichern
    with open(args.out, 'w') as out_vcf:
        out_vcf.write("##fileformat=VCFv4.2\n")
        out_vcf.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n")
        out_vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n")
        out_vcf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
        out_vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        out_vcf.write("##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">\n")
        out_vcf.write("##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another\">\n")
        out_vcf.write("##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">\n")
        out_vcf.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n")
        out_vcf.write("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n")
        out_vcf.write("##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description=\"Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)\">\n")
        out_vcf.write("##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">\n")
        out_vcf.write("##contig=<ID=chr1,length=249250621>\n")
        out_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUnnamedSample\n")
        
        if cn_state < 0.5:
            out_vcf.write(f"chr1\t{POS_Allele}\t.\t{REF_Allele}\t{ALT_Allele}\t255.93\t.\tAC=2;AF=1;AN=2;DP=60;ExcessHet=0;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=30.03;SOR=1.329\tGT:AD:DP:GQ:PL\t1/1:30,30:30:18:270,18,0\n")
        elif cn_state < 1.5:
            out_vcf.write(f"chr1\t{POS_Allele}\t.\t{REF_Allele}\t{ALT_Allele}\t255.93\t.\tAC=2;AF=1;AN=2;DP=60;ExcessHet=0;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=30.03;SOR=1.329\tGT:AD:DP:GQ:PL\t0/1:0,60:60:18:270,18,0\n")
        else:
            out_vcf.write(f"chr1\t{POS_Allele}\t.\t{REF_Allele}\t{ALT_Allele}\t255.93\t.\tAC=2;AF=1;AN=2;DP=60;ExcessHet=0;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=30.03;SOR=1.329\tGT:AD:DP:GQ:PL\t0/0:60,0:60:18:270,18,0\n")


if __name__ == "__main__":
    main()

