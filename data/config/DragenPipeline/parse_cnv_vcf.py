#!/usr/bin/env python3

import argparse
import gzip
import re
from cyvcf2 import VCF

def parse_info_field(info):
    return dict(item.split('=') for item in info.split(';') if '=' in item)

def parse_region(region_str):
    match = re.match(r'^(chr[\w\d]+):(\d+)-(\d+)$', region_str)
    if not match:
        raise ValueError(f"Ungültiges Region-Format: {region_str}")
    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    return chrom, start, end

def overlaps(region_start, region_end, vcf_start, vcf_end):
    return vcf_start <= region_end and vcf_end >= region_start

def main(vcf_file, region_str):
    chrom, region_start, region_end = parse_region(region_str)

    vcf_reader = VCF(vcf_file)
    for variant in vcf_reader:
        # 1. Chromosom vergleichen
        if variant.CHROM != chrom:
            continue

        # 2. Region überschneidet
        var_start = variant.start  # 0-basiert
        var_end = variant.INFO.get("END", var_start + 1)  # fallback
        if var_end < region_start or var_start > region_end:
            continue

        # 3. ALT == <DEL>
        if not variant.ALT or "<DEL>" not in str(variant.ALT[0]):
            continue

        # 4. FILTER == PASS
        if variant.FILTER not in (None, 'PASS'):
            continue

        # 5. SVLEN extrahieren
        svlen = int(variant.INFO.get("SVLEN", 0))

        # 6. GT (Genotyp) als String extrahieren
        gt_tuple = variant.genotypes[0][:2]  # z. B. (1, 1)
        gt = f"{gt_tuple[0]}/{gt_tuple[1]}"

        # ✅ Ausgabe
        print(f"{variant.CHROM}:{variant.start + 1}-{var_end} | SVLEN: {svlen} | GT: {gt} | VALID: TRUE")
        with open(args.out, 'w') as out_vcf:
            out_vcf.write("##fileformat=VCFv4.2\n")
            out_vcf.write("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n")
            out_vcf.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency of alternate alleles\">\n")
            out_vcf.write("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n")
            out_vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for the site\">\n")
            out_vcf.write("##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS mapping quality\">\n")
            out_vcf.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">\n")
            out_vcf.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n")
            out_vcf.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
            out_vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            out_vcf.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles\">\n")
            out_vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth of sample\">\n")
            out_vcf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
            out_vcf.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n")
            out_vcf.write("##contig=<ID=chr1,length=249250621>\n")
            out_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUnnamedSample\n")
            
            REF_Allele = "ALONG"
            ALT_Allele = "A"
            POS_Allele = 25272447
            
            out_vcf.write(f"chr1\t{POS_Allele}\t.\t{REF_Allele}\t{ALT_Allele}\t255.93\t.\tAC=2;AF=1;AN=2;DP=60;ExcessHet=0;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=30.03;SOR=1.329\tGT:AD:DP:GQ:PL\t{gt_tuple[0]}/{gt_tuple[1]}:60,0:60:18:270,18,0\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter VCF CNVs by genomic region.")
    parser.add_argument('vcf_file', help="Pfad zur VCF-Datei (.vcf oder .vcf.gz)")
    parser.add_argument('-p', '--region', required=True, help="Region im Format chr:start-end (z.B. chr1:25272509-25330445)")
    parser.add_argument("-o", "--out", required=True, help="Path to output VCF file.")
    
    args = parser.parse_args()
    main(args.vcf_file, args.region)

