import os
import subprocess
import concurrent.futures

input_path = os.path.abspath("./")
output_path = os.path.abspath("./")
os.makedirs(output_path, exist_ok=True)

# OmegaPlus parameters
minwin_values = [7500]  # could also be: minwin_values = [5000, 7500, 10000], if more than one value is tested
maxwin_values = [75000]  # could also be: maxwin_values = [50000, 75000, 100000], if more than one value is tested

# Chromosome lengths of Arabis nemo ref genome
chromosome_lengths = {
    'chr1': 31339903,
    'chr2': 26594995,
    'chr3': 30118065,
    'chr4': 41679211,
    'chr5': 26526041,
    'chr6': 22121730,
    'chr7': 26822521,
    'chr8': 39891204
}


base_patterns = [
    "combined_arabis_polar_phased_filtered_final_allSNPs_ref_nemo_chr{chr_num}.RUSrem.ArabisNemo.vcf",
    "combined_arabis_polar_phased_filtered_final_allSNPs_ref_nemo_chr{chr_num}.RUSrem.ArabisSagit.vcf"
]

# match the patterns for chromosomes 1 to 8
vcf_files = [
    f for f in os.listdir(input_path) if f.endswith(".vcf") and
    any(f.startswith(pattern.format(chr_num=i).split("chr{chr_num}")[0]) for pattern in base_patterns for i in range(1, 9))
]

# run OmegaPlus
def run_OmegaPlus(vcf_file, grid, minwin, maxwin, outputdir):
    parts = vcf_file.split('_')
    seed = parts[7]  # seed from file name
    label = parts[8].split('.')[0]  # label from file name
    chr_part = [part for part in parts if part.startswith('chr')][0]  # chromosome part from file name
    chr_num = chr_part.split('.')[0]  # extracting chromosome number (e.g., 'chr1')
    model = 'ArabisNemo' if 'ArabisNemo' in vcf_file else 'ArabisSagit'
    vcf_file_path = os.path.join(input_path, vcf_file)  # path to the VCF file
    length = chromosome_lengths[chr_num]  # length of the chromosome
    command = f"/home/lukas/software/omegaplus/OmegaPlus -name {seed}_{label}_{model}_grid{grid}_minwin{minwin}_maxwin{maxwin}_{chr_num} -input {vcf_file_path} -grid {grid} -minwin {minwin} -maxwin {maxwin} -length {length} -seed 2"
    print(f"Running OmegaPlus command: {command}")

    # Run OmegaPlus
    result = subprocess.run(command, cwd=outputdir, shell=True, capture_output=True, text=True)

    # debug
    if result.stdout:
        print(f"OmegaPlus output for {label}:\n{result.stdout}")
    if result.stderr:
        print(f"OmegaPlus error for {label}:\n{result.stderr}")

    # Check / debug
    results_file = os.path.join(outputdir, f"OmegaPlus_Report.{seed}_{label}_{model}_grid{grid}_minwin{minwin}_maxwin{maxwin}_{chr_num}")
    if not os.path.exists(results_file):
        raise FileNotFoundError(f"Expected OmegaPlus results file not found: {results_file}")

# processing each VCF file
def process_vcf_file(vcf_file):
    try:
        parts = vcf_file.split('_')
        chr_part = [part for part in parts if part.startswith('chr')][0]  # get chromosome part from file name
        chr_num = chr_part.split('.')[0]  
        chromosome_length = chromosome_lengths[chr_num]
        grid = round(chromosome_length / 100000)  # calculate grid value based on chromosome length

        for minwin in minwin_values:
            for maxwin in maxwin_values:
                run_OmegaPlus(vcf_file, grid, minwin, maxwin, output_path)
    except Exception as e:
        print(f"Error processing {vcf_file}: {e}")

# Number of threads
num_threads = 2 

# run in parallel
with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
    futures = {executor.submit(process_vcf_file, vcf_file): vcf_file for vcf_file in vcf_files}
    for future in concurrent.futures.as_completed(futures):
        future.result()

print("Processing complete. OmegaPlus runs complete.")
