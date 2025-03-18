import os
import subprocess
import sys
# Define Directories and Paths
NGS_DIR = "/home/user/Desktop/NGS_Final"
REFERENCE_DIR = os.path.join(NGS_DIR, "reference_genome")
REFERENCE_GENOME = os.path.join(REFERENCE_DIR, "hg38.fa")
KNOWN_SITES = os.path.join(REFERENCE_DIR, "Homo_sapiens_assembly38.known_indels.vcf.gz")
TRIMMOMATIC_JAR = os.path.join(NGS_DIR, "Trimmomatic-0.39/trimmomatic-0.39.jar")
ADAPTERS = os.path.join(NGS_DIR, "Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa")
GATK = "/usr/local/bin/gatk/gatk"
PICARD = "java -jar /home/user/Downloads/picard.jar"
VARIANT_FILTER_SCRIPT = "/home/user/Desktop/NGS_Final/variant_filteration.py"
THREADS = 50  # Adjust as needed for system performance

# Required directories
REQUIRED_DIRS = ["qc_reports", "trimmed_outputs", "aligned", "results"]

def execute_command(command, step_desc=""):
    """Run a shell command and handle errors."""
    try:
        print(f"🔹 Running: {step_desc}...")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error in {step_desc}: {command}")
        print(f"⚠️ {e}")
        sys.exit(1)  # Exit script on failure

def ensure_reference_index():
    """Ensure reference genome is indexed correctly (FASTA index & dictionary)."""
    fai_file = f"{REFERENCE_GENOME}.fai"
    dict_file = f"{REFERENCE_GENOME}.dict"

    if not os.path.exists(fai_file):
        print("⚠️ FASTA index (.fai) missing. Creating index...")
        execute_command(f"samtools faidx {REFERENCE_GENOME}", "Generating FASTA Index")

    if not os.path.exists(dict_file):
        print("⚠️ FASTA dictionary (.dict) missing. Creating dictionary...")
        execute_command(f"{GATK} CreateSequenceDictionary -R {REFERENCE_GENOME} -O {dict_file}", "Generating FASTA Dictionary")

def ensure_bwa_index():
    """Ensure BWA index files exist for the reference genome."""
    required_index_files = [
        f"{REFERENCE_GENOME}.amb",
        f"{REFERENCE_GENOME}.ann",
        f"{REFERENCE_GENOME}.bwt",
        f"{REFERENCE_GENOME}.pac",
        f"{REFERENCE_GENOME}.sa",
    ]
    if all(os.path.exists(f) for f in required_index_files):
        print("✅ BWA index files found.")
    else:
        print("⚠️ BWA index files missing. Creating index...")
        execute_command(f"bwa index {REFERENCE_GENOME}", "BWA Indexing")

def add_read_groups(patient, dedup_bam):
    """Check if the BAM file has read groups, if not, add them."""
    print(f"🔍 Checking read groups for {patient}...")
    check_rg_cmd = f"samtools view -H {dedup_bam} | grep '@RG'"
    result = subprocess.run(check_rg_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if not result.stdout.strip():
        print(f"⚠️ No read groups found in {dedup_bam}. Adding read groups...")
        dedup_with_rg_bam = dedup_bam.replace(".bam", "_with_rg.bam")
        add_rg_cmd = (
            f"{PICARD} AddOrReplaceReadGroups "
            f"-I {dedup_bam} "
            f"-O {dedup_with_rg_bam} "
            f"-RGID 1 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM {patient}"
        )
        execute_command(add_rg_cmd, "Adding Read Groups")
        return dedup_with_rg_bam
    else:
        print(f"✅ Read groups already present in {dedup_bam}.")
        return dedup_bam

def process_patient(patient):
    """Run variant calling pipeline for a patient."""
    print(f"\n🚀 Processing Patient: {patient}")

    # Ensure Reference Indexes
    ensure_reference_index()
    ensure_bwa_index()

    # Directories
    patient_dir = os.path.join(NGS_DIR, patient)
    qc_output = os.path.join(patient_dir, "qc_reports")
    trim_output = os.path.join(patient_dir, "trimmed_outputs")
    align_output = os.path.join(patient_dir, "aligned")
    results_output = os.path.join(patient_dir, "results")

    # Ensure directories exist
    os.makedirs(qc_output, exist_ok=True)
    os.makedirs(trim_output, exist_ok=True)
    os.makedirs(align_output, exist_ok=True)
    os.makedirs(results_output, exist_ok=True)

    # Files
    r1_file = os.path.join(patient_dir, f"{patient}_R1.fastq.gz")
    r2_file = os.path.join(patient_dir, f"{patient}_R2.fastq.gz")
    trimmed_r1 = os.path.join(trim_output, f"{patient}_R1_paired.fastq.gz")
    trimmed_r2 = os.path.join(trim_output, f"{patient}_R2_paired.fastq.gz")

    # ✅ Step 1: Quality Control & Trimming
    execute_command(f"fastqc -t {THREADS} -o {qc_output} {r1_file} {r2_file}", "FastQC Quality Control")
    execute_command(
        f"java -jar {TRIMMOMATIC_JAR} PE -threads {THREADS} {r1_file} {r2_file} "
        f"{trimmed_r1} {trim_output}/{patient}_R1_unpaired.fastq.gz "
        f"{trimmed_r2} {trim_output}/{patient}_R2_unpaired.fastq.gz "
        f"ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
        "Adapter Trimming & Quality Filtering"
    )

    # ✅ Step 2: Read Alignment
    sam_file = os.path.join(align_output, f"{patient}.sam")
    execute_command(f"bwa mem -t {THREADS} {REFERENCE_GENOME} {trimmed_r1} {trimmed_r2} > {sam_file}", "BWA Alignment")

    # ✅ Step 3: Processing BAM Files
    dedup_bam = os.path.join(align_output, f"{patient}_dedup.bam")
    recalibrated_bam = os.path.join(align_output, f"{patient}_recalibrated.bam")
    raw_vcf = os.path.join(results_output, "raw_variants.vcf")
    bamout_file = os.path.join(results_output, f"{patient}_bamout.bam")

    execute_command(f"samtools sort {sam_file} -o {dedup_bam}", "Sorting BAM")
    execute_command(f"samtools index {dedup_bam}", "Indexing BAM")
    
    dedup_with_rg_bam = add_read_groups(patient, dedup_bam)

    execute_command(f"{GATK} BaseRecalibrator -R {REFERENCE_GENOME} -I {dedup_with_rg_bam} --known-sites {KNOWN_SITES} -O {results_output}/recal_data.table", "Base Recalibration")
    execute_command(f"{GATK} ApplyBQSR -R {REFERENCE_GENOME} -I {dedup_with_rg_bam} --bqsr-recal-file {results_output}/recal_data.table -O {recalibrated_bam}", "Applying BQSR")
    execute_command(
        f"{GATK} HaplotypeCaller "
        f"--java-options '-Xmx16g' "
        f"-R {REFERENCE_GENOME} "
        f"-I {recalibrated_bam} "
        f"-O {raw_vcf} "
        f"-ERC GVCF "
        f"-bamout {bamout_file} "
        f"--native-pair-hmm-threads {THREADS} "
        f"-A StrandBiasBySample "
        f"-A FisherStrand "
        f"-A StrandOddsRatio "
        f"-A GenotypeSummaries "
        f"-A RMSMappingQuality "
        f"-A MappingQualityRankSumTest "
        f"-A ReadPosRankSumTest "
        f"-A QualByDepth "
        f"-A DepthPerAlleleBySample "
        f"-A DepthPerSampleHC "
        f"-A Coverage "
        f"-A SampleList "
        f"-G StandardAnnotation "
        f"-G AS_StandardAnnotation "
        f"-G StandardHCAnnotation"
    )
    print(f"✅ Processing Completed for {patient}")

    # ✅ Step 4: Run Variant Filtering Script
    print(f"🔹 Running Variant Filtering for {patient}...")
    execute_command(f"python3 {VARIANT_FILTER_SCRIPT} {patient}", "Running Variant Filtering")

if __name__ == "__main__":
    process_patient(sys.argv[1])

