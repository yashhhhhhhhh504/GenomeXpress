import os
import subprocess
import sys
import pandas as pd

# Define Directories and Paths
NGS_DIR = "/home/user/Desktop/NGS_Final"
REFERENCE_DIR = os.path.join(NGS_DIR, "reference_genome")
REFERENCE_GENOME = os.path.join(REFERENCE_DIR, "hg38.fa")
ANNOVAR_DIR = os.path.join(NGS_DIR, "annovar")
HUMANDB_DIR = os.path.join(ANNOVAR_DIR, "humandb")
TABLE_ANNOVAR = os.path.join(ANNOVAR_DIR, "table_annovar.pl")
CONVERT_ANNOVAR = os.path.join(ANNOVAR_DIR, "convert2annovar.pl")

# List of databases for annotation
DATABASES = [
    "refGene", "avsnp150", "gnomad_genome", "esp6500siv2_all",
    "clinvar_20210501", "revel", "dbnsfp42a", "dbscsnv11",
    "genomicSuperDups", "gwasCatalog", "kgXref"
]

def execute_command(command, description):
    """Execute a shell command with error handling."""
    try:
        print(f"\n🔹 {description}...")
        subprocess.run(command, shell=True, check=True)
        print(f"✅ {description} completed.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error during {description}: {e}")
        sys.exit(1)

def annotate_variants(patient_name):
    """Annotate filtered VCF variants using ANNOVAR."""
    results_dir = os.path.join(NGS_DIR, patient_name, "results")
    filtered_vcf = os.path.join(results_dir, "filtered_variants.vcf.gz")

    if not os.path.exists(filtered_vcf):
        print(f"❌ Error: Filtered VCF file not found at {filtered_vcf}")
        sys.exit(1)

    # Convert VCF to ANNOVAR input format
    avinput_file = os.path.join(results_dir, "filtered_variants.avinput")
    execute_command(
        f"perl {CONVERT_ANNOVAR} -format vcf4 {filtered_vcf} > {avinput_file}",
        "Converting VCF to ANNOVAR input format"
    )

    # Run ANNOVAR annotation
    output_prefix = os.path.join(results_dir, "annotated_variants")
    annotation_command = (
        f"perl {TABLE_ANNOVAR} {avinput_file} {HUMANDB_DIR} "
        f"-buildver hg38 -out {output_prefix} -remove "
        f"-protocol {','.join(DATABASES)} "
        f"-operation {','.join(['g' if db == 'refGene' else 'f' for db in DATABASES])} "
        f"-nastring . --csvout"
    )
    execute_command(annotation_command, "Running ANNOVAR Annotation with CSV output")

    # Convert annotated results into CSV
    annovar_output_file = output_prefix + ".hg38_multianno.txt"
    csv_output_file = output_prefix + ".csv"
    convert_annovar_output_to_csv(annovar_output_file, csv_output_file, filtered_vcf)

    print(f"\n🎉 Annotation completed successfully for {patient_name}!")

def convert_annovar_output_to_csv(annovar_output_file, csv_output_file, vcf_file):
    """Convert ANNOVAR output and include VCF fields into CSV format."""
    print("\n🔹 Converting ANNOVAR output to CSV format...")

    # Read ANNOVAR annotation results
    annovar_df = pd.read_csv(annovar_output_file, sep="\t", header=0, low_memory=False)

    # Extract VCF information (DP, AD, AF, GT) and merge with ANNOVAR results
    vcf_data = extract_vcf_info(vcf_file)
    final_df = pd.merge(annovar_df, vcf_data, on=["Chr", "Start", "Ref", "Alt"], how="left")

    # Save final merged output
    final_df.to_csv(csv_output_file, index=False)
    print(f"✅ ANNOVAR output saved as CSV: {csv_output_file}")

def extract_vcf_info(vcf_file):
    """Extract key fields from the filtered VCF to include in the final CSV output."""
    vcf_data = []
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip headers
            fields = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, info, _, genotype = fields[:10]

            # Extract key INFO fields (DP, AD, AF, etc.)
            info_dict = {key: value for key, value in (entry.split("=") for entry in info.split(";") if "=" in entry)}
            depth = info_dict.get("DP", ".")
            allele_freq = info_dict.get("AF", ".")
            allele_depth = info_dict.get("AD", ".")

            # Extract genotype information
            genotype_info = genotype.split(":")
            gt = genotype_info[0] if len(genotype_info) > 0 else "."

            vcf_data.append([chrom, pos, ref, alt, depth, allele_freq, allele_depth, gt])

    return pd.DataFrame(vcf_data, columns=["Chr", "Start", "Ref", "Alt", "Depth", "Allele_Freq", "Allele_Depth", "Genotype"])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 annotation.py <patient_name>")
        sys.exit(1)

    patient_name = sys.argv[1]
    annotate_variants(patient_name)
