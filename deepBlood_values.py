#!/bin/python3
import json
import sys

# JSON-Dateien aus dem Muster --json lesen
json_files = sys.argv[1:]

# Ausgabe der Header
print("Sample_ID\tLocus\tCall_0_Names\tCall_1_Names\tCall_0_Score\tCall_0_weak_Score\tHaplotypes_0_Variations\tHaplotypes_1_Variations\tno. of relevant_variants\tno. of coverage_failed_variants\tcoverage_failed_variants\trequired_coverage\tmean_coverage_CDS\tmean_coverage_exons")

# Iterieren über die JSON-Dateien
for json_file in json_files:
    # Laden der JSON-Daten
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Iterieren über die Loci
    for locus, locus_data in data["loci"].items():
        sample_id = data.get("sample_id", "")
        coverage_threshold = data["parameters"].get("--coverage", "")
        
        if "calls" in locus_data and len(locus_data["calls"]) > 0:
            for call_data in locus_data["calls"]:
                call_0_names = ",".join(call_data["alleles"][0]["names"])
                call_0_score = "{:.3f}".format(call_data["score"])
                call_0_weak_score = "{:.3f}".format(call_data["weak_score"])
                
                if len(call_data["alleles"]) > 1:
                    call_1_names = ",".join(call_data["alleles"][1]["names"])
                else:
                    call_1_names = "-"
                
                # Haplotypes extrahieren
                if "haplotypes" in call_data:
                    if len(call_data["haplotypes"]["genotypes"]) > 0:
                        haplotypes_0_genotypes = call_data["haplotypes"]["genotypes"][0]
                        if isinstance(call_data["haplotypes"]["genotypes"][0], list):
                            haplotypes_0_variations = ",".join([haplotype["variation"] for haplotype in haplotypes_0_genotypes if haplotype and isinstance(haplotype, dict)])
                        else:
                            haplotypes_0_variations = ""
                    else:
                        haplotypes_0_variations = ""
                    if len(call_data["haplotypes"]["genotypes"]) > 1:
                        haplotypes_1_genotypes = call_data["haplotypes"]["genotypes"][1]
                        if isinstance(call_data["haplotypes"]["genotypes"][1], list):
                            haplotypes_1_variations = ",".join([haplotype["variation"] for haplotype in haplotypes_1_genotypes if haplotype and isinstance(haplotype, dict)])
                        else:
                            haplotypes_1_variations = ""
                    else:
                        haplotypes_1_variations = ""
                else:
                    haplotypes_0_variations = ""
                    haplotypes_1_variations = ""
    
                # Coverage failed variants extrahieren
                coverage_failed_variants = "-"
                no_coverage_failed_variants = "0"
                no_relevant_variants = "0"
                if "relevant_variations" in locus_data and isinstance(locus_data["relevant_variations"], list):
                    no_relevant_variants = len(locus_data["relevant_variations"])
                if "coverage_failed_variants" in locus_data and isinstance(locus_data["coverage_failed_variants"], list):
                    coverage_failed_variants = ",".join(locus_data["coverage_failed_variants"])
                    no_coverage_failed_variants = len(locus_data["coverage_failed_variants"])
                
                mean_coverage_cds = "-"
                mean_coverage_exons = "-"
                if "mean_coverage" in locus_data and locus_data["mean_coverage"]:
                    if "CDS" in locus_data["mean_coverage"]:
                        mean_coverage_cds = "{:.1f}".format(str(locus_data["mean_coverage"]["CDS"]))
                    if "exons" in locus_data["mean_coverage"]:
                        mean_coverage_exons = ",".join([f"{round(exon, 1)}" for exon in locus_data["mean_coverage"]["exons"]])

                # Ausgabe der Zeile
                print(f"{sample_id}\t{locus}\t{call_0_names}\t{call_1_names}\t{call_0_score}\t{call_0_weak_score}\t{haplotypes_0_variations}\t{haplotypes_1_variations}\t{no_relevant_variants}\t{no_coverage_failed_variants}\t{coverage_failed_variants}\t{coverage_threshold}\t{mean_coverage_cds}\t{mean_coverage_exons}")
        else:
            # Wenn keine Anrufe vorhanden sind, geben Sie leere Werte aus
            print(f"{sample_id}\t{locus}\t\t\t\t\t\t\t\t\t")

