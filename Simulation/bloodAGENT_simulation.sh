#!/bin/bash

ROOTPATH="./"
OUTPATH="/tmp"
dropout_freq=5
crack_freq=20


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ROOTPATH}/external/htslib:${ROOTPATH}/external/libBigWig


# rm ${OUTPATH}/result_${dropout_freq}_percent_dropouts_${crack_freq}_percent_haplotype_cracks.txt

# Datei mit den Daten
allele_table="${ROOTPATH}/data/config/genotype_to_phenotype_annotation_TGS.dat"

# Temporäre Dateien
unique_file=$(mktemp)
filtered_file=$(mktemp)
vcf_file=$(mktemp --suffix=".vcf")
json_file=$(mktemp --suffix=".json")


for i in {1..16384}; do


# 1. Spalte 2 einzigartig machen und zufälligen Eintrag auswählen
sed 1D "$allele_table" | grep -v -P "^#" | cut -f2  | sort | uniq | shuf -n 1 > "$unique_file"
selected_key=$(cat "$unique_file")
####################################
#selected_key="JK"

# 2. Zeilen mit demselben Wert in Spalte 2 filtern
grep -v -P "^#" "$allele_table" | awk -v key="$selected_key" '$2 == key'  > "$filtered_file"

# 3. Zufällige zwei Einträge aus den gefilterten Zeilen auswählen (mit Wiederholungen)
selected_entries=$(shuf -n 2 "$filtered_file" | cut -f 4 | tr '\n' ' ')

allele_one=$(echo $selected_entries | cut -f 1 -d ' ')
allele_two=$(echo $selected_entries | cut -f 2 -d ' ')

# Ausgabe der Ergebnisse
cat ${ROOTPATH}/data/simulation_vcf_header.vcf > $vcf_file

${ROOTPATH}/dist/Debug/GNU-Linux/deepblood --job vcf \
          --variants ${ROOTPATH}/data/config/variation_annotation_TGS.dat \
          --gt2pt ${ROOTPATH}/data/config/genotype_to_phenotype_annotation_TGS.dat -a "$allele_one" -b "$allele_two" \
          --phased \
          --dropout ${dropout_freq} --crack ${crack_freq} 1>> $vcf_file


${ROOTPATH}/dist/Debug/GNU-Linux/deepblood -j phenotype \
                  --target ${ROOTPATH}/data/config/exonic_annotation.hg38.BGStarget.txt \
                  --variants ${ROOTPATH}/data/config/variation_annotation_TGS.dat \
                  --gt2pt ${ROOTPATH}/data/config/genotype_to_phenotype_annotation_TGS.dat \
                  --vcf $vcf_file \
                  --bigwig ${ROOTPATH}/data/config/fake_50x_coverage_at_target.hg38.bw \
                  --coverage 0 \
                  --verbose 2 \
                  --scoreRange 0.99 \
                  --out $json_file \
                  --build hg38 --id "Simulation" --locus "$selected_key"
                  
RESULT=$(python ${ROOTPATH}/deepBlood_values.py $json_file | grep Simulation | cut -f 2-4)

echo "$selected_key $selected_entries $RESULT"
echo "$selected_key $selected_entries $RESULT" | awk '{
    # Sortiere die Spalten 2 und 3
    if ($2 < $3) {
        sorted1 = $2 " " $3
    } else {
        sorted1 = $3 " " $2
    }
    
    # Sortiere die Spalten 5 und 6
    if ($5 < $6) {
        sorted2 = $5 " " $6
    } else {
        sorted2 = $6 " " $5
    }
    
    # Drucke die Zeile mit sortierten Spalten
    printf "%s\t%s\t%s\t%s\n", $1, sorted1, $4, sorted2
}' >> ${OUTPATH}/result_${dropout_freq}_percent_dropouts_${crack_freq}_percent_haplotype_cracks.txt

done

# Aufräumen temporärer Dateien
rm "$unique_file" "$filtered_file" "$vcf_file" "$json_file"


