import os
import subprocess
import sys
import multiprocessing
NGS_DIR = "/home/user/Desktop/NGS_Final"
REFERENCE_GENOME = os.path.join(NGS_DIR, "reference_genome", "hg38.fa")
GATK = "/usr/local/bin/gatk/gatk"
# **Detect available CPU threads & assign an optimal number**
TOTAL_THREADS = os.cpu_count() if os.cpu_count() else 8  # Get CPU count, fallback to 8
THREADS = max(4, int(TOTAL_THREADS * 0.75))  # Use 75% of available cores, minimum 4

def execute_command(command, description):
    """Run a shell command and handle errors."""
    try:
        print(f"\n🔹 {description}...")
        subprocess.run(command, shell=True, check=True)
        print(f"✅ {description} completed.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error during {description}: {e}")
        sys.exit(1)

def validate_vcf(vcf_path):
    """Check if a VCF file is valid and properly formatted."""
    if not os.path.exists(vcf_path):
        print(f"❌ Error: VCF file not found at {vcf_path}")
        sys.exit(1)

    # Check if the VCF file is compressed properly
    check_command = f"bcftools view {vcf_path} > /dev/null 2>&1"
    if subprocess.call(check_command, shell=True) != 0:
        print(f"⚠️ Warning: {vcf_path} is not a valid VCF. Recompressing...")
        recompress_command = f"bgzip -c {vcf_path} > {vcf_path}.gz && tabix -p vcf {vcf_path}.gz"
        execute_command(recompress_command, f"Recompressing and indexing {vcf_path}")
        return f"{vcf_path}.gz"
    
    return vcf_path

def process_patient(patient):
    results_dir = os.path.join(NGS_DIR, patient, "results")
    raw_vcf = os.path.join(results_dir, "raw_variants.vcf.gz")
    combined_vcf = os.path.join(results_dir, "combined_variants.g.vcf.gz")
    genotyped_vcf = os.path.join(results_dir, "genotyped_variants.vcf.gz")
    normalized_vcf = os.path.join(results_dir, "normalized_variants.vcf.gz")
    filtered_vcf = os.path.join(results_dir, "filtered_variants.vcf.gz")
    if not os.path.exists(raw_vcf):
        print(f"❌ Error: Raw variants file not found for {patient} at {raw_vcf}")
        sys.exit(1)
    raw_vcf = validate_vcf(raw_vcf)
    execute_command(
        f"{GATK} --java-options '-Xmx16g' CombineGVCFs "
        f"-R {REFERENCE_GENOME} --variant {raw_vcf} -O {combined_vcf} "
        f"--tmp-dir {results_dir}",
        "Combining GVCFs"
    )
    execute_command(
        f"{GATK} --java-options '-Xmx16g' GenotypeGVCFs "
        f"-R {REFERENCE_GENOME} -V {combined_vcf} -O {genotyped_vcf} "
        f"--native-pair-hmm-threads {THREADS}",
        "Genotyping GVCFs"
    )
    genotyped_vcf = validate_vcf(genotyped_vcf)
    execute_command(
        f"bcftools norm -f {REFERENCE_GENOME} --multiallelics -any "
        f"-O z -o {normalized_vcf} {genotyped_vcf}",
        "Normalizing Variants"
    )
    execute_command(
        f"tabix -p vcf {normalized_vcf}",
        "Indexing Normalized VCF"
    )
    execute_command(
        f"{GATK} VariantFiltration -R {REFERENCE_GENOME} -V {normalized_vcf} -O {filtered_vcf} "
        f"--filter-expression 'QD < 2.0' --filter-name 'QD2' "
        f"--filter-expression 'QUAL < 30.0' --filter-name 'QUAL30' "
        f"--filter-expression 'SOR > 3.0' --filter-name 'SOR3' "
        f"--filter-expression 'FS > 60.0' --filter-name 'FS60' "
        f"--filter-expression 'MQ < 40.0' --filter-name 'MQ40' "
        f"--filter-expression 'DP < 10' --filter-name 'DP10'",
        "Applying Variant Filtration"
    )
    print(f"\n🎉 Variant Processing Completed for {patient}")
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 variant_filteration.py <patient_name>")
        sys.exit(1)
    patient_name = sys.argv[1]
    process_patient(patient_name)