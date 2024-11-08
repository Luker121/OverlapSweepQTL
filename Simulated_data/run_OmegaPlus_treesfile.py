import msprime
import tskit
import pyslim
import numpy as np
import os
import concurrent.futures
import warnings
import subprocess
import time
from msprime import TimeUnitsMismatchWarning


warnings.simplefilter('ignore', TimeUnitsMismatchWarning)

# paths
input_path = os.path.abspath("./")
output_path = os.path.abspath("./OmegaPlus/")
os.makedirs(output_path, exist_ok=True)
tree_files = [f for f in os.listdir(input_path) if f.endswith(".trees")]

# Function to read a parameter file from SLiM simulations and extract the NAnc and mAncNS values
def get_parameters(param_file):
    NAnc = None
    mAncNS = None
    try:
        with open(param_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("NAnc="):
                    NAnc = int(line.split('=')[1].strip())
                elif line.startswith("mAncNS="):
                    mAncNS = float(line.split('=')[1].strip())
        if NAnc is not None and mAncNS is not None:
            return NAnc, mAncNS
    except UnicodeDecodeError:
        with open(param_file, 'r', encoding='latin1') as f:
            for line in f:
                if line.startswith("NAnc="):
                    NAnc = int(line.split('=')[1].strip())
                elif line.startswith("mAncNS="):
                    mAncNS = float(line.split('=')[1].strip())
        if NAnc is not None and mAncNS is not None:
            return NAnc, mAncNS
    raise ValueError(f"Parameters NAnc or mAncNS not found in {param_file}")

# Function to get individuals from a specific population
def get_individuals_from_population(ts, alive, population):
    return [ind for ind in alive if ts.node(ts.individual(ind).nodes[0]).population == population]

# Function to get nodes from individuals
def get_nodes_from_individuals(ts, individuals):
    nodes = []
    for ind in individuals:
        nodes.extend(ts.individual(ind).nodes)
    return nodes

# Function to run OmegaPlus
def run_OmegaPlus(ms_file_path, grid, minwin, maxwin, length, outputdir, seed, label):
    command = f"/home/ubuntu/omegaplus/OmegaPlus -name {seed}_{label}_grid{grid}_minwin{minwin}_maxwin{maxwin} -input {ms_file_path} -grid {grid} -minwin {minwin} -maxwin {maxwin} -length {length}"
    print(f"Running OmegaPlus command: {command}")
    
    # Run OmegaPlus and capture output and error messages
    result = subprocess.run(command, cwd=outputdir, shell=True, capture_output=True, text=True)
    
    # Print output and error messages
    if result.stdout:
        print(f"OmegaPlus output for {label}:\n{result.stdout}")
    if result.stderr:
        print(f"OmegaPlus error for {label}:\n{result.stderr}")
    
    # Check if the results file was created
    results_file = os.path.join(outputdir, f"OmegaPlus_Report.{seed}_{label}_grid{grid}_minwin{minwin}_maxwin{maxwin}")
    if not os.path.exists(results_file):
        raise FileNotFoundError(f"Expected OmegaPlus results file not found: {results_file}")
    
    # Delete the .ms file to save space
    os.remove(ms_file_path)

# Function to extract seed from file name
def extract_seed(file_name):
    parts = file_name.split('_')
    for i, part in enumerate(parts):
        if part == 'seed':
            return parts[i+1]
    return None

# Function to process each tree file
def process_tree_file(tree_file, file_count, total_files):
    try:
        # Load the tree sequence
        orig_ts = tskit.load(os.path.join(input_path, tree_file))
        print(f"Processing {tree_file}: {orig_ts.num_trees} trees, {orig_ts.num_mutations} mutations.")
        
        # Get the corresponding parameter file
        param_file = tree_file.replace("overlay.trees", "parameters.txt")
        param_file_path = os.path.join(input_path, param_file)
        
        # Check if the parameter file exists
        if not os.path.exists(param_file_path):
            print(f"Parameter file {param_file_path} not found, skipping.")
            return

        # Get the ancestral Ne and mAncNS values from the parameter file
        try:
            ancestral_Ne, _ = get_parameters(param_file_path)
        except ValueError as e:
            print(e)
            print(f"Skipping {tree_file} due to missing parameters.")
            return

        # Recapitate and add mutations
        recap_ts = pyslim.recapitate(orig_ts, recombination_rate=5.56e-7, ancestral_Ne=ancestral_Ne)
        ts = msprime.sim_mutations(
            recap_ts,
            rate=7e-8,
            model=msprime.SLiMMutationModel(type=0),
            keep=False,
            discrete_genome=False
        )
        
        print(f"The tree sequence now has {ts.num_trees} trees,\n"
              f"and {ts.num_mutations} mutations.")

        # Get individuals alive at the end of the simulation
        alive = pyslim.individuals_alive_at(recap_ts, 0)
        print(f"There are {len(alive)} individuals alive in the final generation.")

        # Sample individuals from specified populations
        rng = np.random.default_rng(seed=1)
        populations_to_sample = {
            "sagit": {5: 15, 4: 10},  # sagit sym and allo
            "nemo": {7: 6, 6: 6}  # nemo sym and allo
        }

        combined_samples = {
            "sagit": [],
            "nemo": []
        }

        for label, pops in populations_to_sample.items():
            for pop, num_samples in pops.items():
                pop_individuals = get_individuals_from_population(ts, alive, pop)
                if len(pop_individuals) >= num_samples:
                    sampled = rng.choice(pop_individuals, num_samples, replace=False)
                    combined_samples[label].extend(sampled)
                else:
                    print(f"Not enough individuals in population {pop} to sample {num_samples} individuals.")
                    combined_samples[label].extend(pop_individuals)

        # Print sampled individuals for debug
        print(f"Combined samples: {combined_samples}")

        # Extract seed from the tree file name
        seed = extract_seed(tree_file)

        # Create .ms files for combined samples
        grid = 40
        minwin = 5000
        maxwin = 50000
        length = int(ts.sequence_length)

        for label, individuals in combined_samples.items():
            nodes = get_nodes_from_individuals(ts, individuals)
            simplified_ts = ts.simplify(nodes, filter_sites=False)
            ms_file_path = os.path.join(output_path, f"out_{seed}_{label}.ms")
            with open(ms_file_path, 'w') as ms_file:
                tskit.write_ms(simplified_ts, ms_file)
            print(f"Created .ms file for {label}: {ms_file_path}")
            run_OmegaPlus(ms_file_path, grid, minwin, maxwin, length, output_path, seed, label)
        
        print(f"Processed {file_count}/{total_files} files.")
        
    except Exception as e:
        print(f"Error processing {tree_file}: {e}")

# Number of threads
num_threads = 10 

# Run in parallel
total_files = len(tree_files)
with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
    futures = {executor.submit(process_tree_file, tree_file, i+1, total_files): tree_file for i, tree_file in enumerate(tree_files)}
    for future in concurrent.futures.as_completed(futures):
        future.result()

print("Processing complete. OmegaPlus runs complete.")
