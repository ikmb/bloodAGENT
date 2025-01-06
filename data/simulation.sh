#!/bin/bash

export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/home/mwittig/coding/cpp/MyTools/dist/Debug/GNU-Linux/:/home/mwittig/coding/fremd/htslib:/home/mwittig/coding/fremd/libBigWig

# Datei mit den Daten
allele_table="/home/mwittig/coding/cpp/deepBlood/data/config/genotype_to_phenotype_annotation_TGS.dat"


while true; do
# Temporäre Dateien
unique_file=$(mktemp)
filtered_file=$(mktemp)
vcf_file=$(mktemp --suffix=".vcf")
json_file=$(mktemp --suffix=".json")



# 1. Spalte 2 einzigartig machen und zufälligen Eintrag auswählen
cut -f2 "$allele_table" | sort | uniq | shuf -n 1 > "$unique_file"
selected_key=$(cat "$unique_file")

# 2. Zeilen mit demselben Wert in Spalte 2 filtern
awk -v key="$selected_key" '$2 == key' "$allele_table" > "$filtered_file"

# 3. Zufällige zwei Einträge aus den gefilterten Zeilen auswählen (mit Wiederholungen)
selected_entries=$(shuf -n 2 "$filtered_file" | cut -f 4 | tr '\n' ' ')

allele_one=$(echo $selected_entries | cut -f 1 -d ' ')
allele_two=$(echo $selected_entries | cut -f 2 -d ' ')

# Ausgabe der Ergebnisse
cat /home/mwittig/coding/cpp/deepBlood/data/simulation_vcf_header.vcf > $vcf_file

/home/mwittig/coding/cpp/deepBlood/dist/Debug/GNU-Linux/deepblood --job vcf \
          --variants /home/mwittig/coding/cpp/deepBlood/data/config/variation_annotation_TGS.dat \
          --gt2pt /home/mwittig/coding/cpp/deepBlood/data/config/genotype_to_phenotype_annotation_TGS.dat -a "$allele_one" -b "$allele_two" \
          --phased \
          --dropout 0 1>> $vcf_file
          

/home/mwittig/coding/cpp/deepBlood/dist/Debug/GNU-Linux/deepblood -j phenotype \
                  --target /home/mwittig/coding/cpp/deepBlood/data/config/exonic_annotation.hg38.BGStarget.txt \
                  --variants /home/mwittig/coding/cpp/deepBlood/data/config/variation_annotation_TGS.dat \
                  --gt2pt /home/mwittig/coding/cpp/deepBlood/data/config/genotype_to_phenotype_annotation_TGS.dat \
                  --vcf $vcf_file \
                  --bigwig /home/mwittig/coding/cpp/deepBlood/data/config/fake_50x_coverage_at_target.hg38.bw \
                  --coverage 0 \
                  --verbose 2 \
                  --scoreRange 0.99 \
                  --insilicovcf \
                  --out $json_file \
                  --build hg38 -k --id "Simulation" --locus "$selected_key"
                  
RESULT=$(deepBlood_values.py $json_file | grep Simulation | cut -f 2-4)

echo "$selected_key $selected_entries $RESULT"
echo "$selected_key $selected_entries $RESULT" >> /home/mwittig/ramDisk/result.txt

done


# Aufräumen temporärer Dateien
rm "$unique_file" "$filtered_file" "$vcf_file" "$json_file"


