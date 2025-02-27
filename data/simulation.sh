#!/bin/bash

#####################################################################
#### Excel help
#### =MAX(SUM(NOT(ISERR(FIND(B1,E1))),NOT(ISERR(FIND(C1,F1)))),SUM(NOT(ISERR(FIND(B1,F1))),NOT(ISERR(FIND(C1,E1)))))
#### =(LEN(LOWER(E1)) - LEN(LOWER(SUBSTITUTE(E1, ","&A1, ""))))/(LEN(A1)+1)+(LEN(LOWER(F1)) - LEN(LOWER(SUBSTITUTE(F1, ","&A1, ""))))/(LEN(A1)+1)
#####################################################################

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work_beegfs/${USER}/software/MyTools/dist/Debug/GNU-Linux/:/work_beegfs/${USER}/software/htslib:/work_beegfs/${USER}/software/libBigWig
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mwittig/coding/cpp/MyTools/dist/Debug/GNU-Linux/:/home/mwittig/coding/fremd/htslib:/home/mwittig/coding/fremd/libBigWig

ROOTPATH="/home/mwittig/coding/cpp/deepBlood"
OUTPATH="/home/mwittig/ramDisk"

# Datei mit den Daten
allele_table="${ROOTPATH}/data/config/genotype_to_phenotype_annotation_TGS.dat"


for drop in 0 5 10 15 20 25 30 35 40 45 50;do

for i in {1..5000}; do
# Temporäre Dateien
unique_file=$(mktemp)
filtered_file=$(mktemp)
vcf_file=$(mktemp --suffix=".vcf")
json_file=$(mktemp --suffix=".json")



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
          --dropout ${drop} --crack 50 1>> $vcf_file


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
                  
RESULT=$(deepBlood_values.py $json_file | grep Simulation | cut -f 2-4)

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
}' >> ${OUTPATH}/result_${drop}_percent_dropout_rate.txt

done
done

# Aufräumen temporärer Dateien
rm "$unique_file" "$filtered_file" "$vcf_file" "$json_file"


