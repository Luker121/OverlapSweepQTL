import pandas as pd
import glob

def read_OmegaPlus_files(directory, models=['ArabisNemo', 'ArabisSagit']):
    df_dict = {}
    for model in models:
        pattern = f"{directory}/OmegaPlus_Report.*_{model}_grid*_minwin*_maxwin*_*"
        files = glob.glob(pattern)
        print(f"Files for model '{model}' found: {len(files)}")
        
        # store dataframes based on unique grid, minwin, and maxwin
        model_dfs = {}
        
        for file in files:
            # Extract grid, minwin, maxwin, and chromosome from the file name
            filename = file.split('/')[-1]  
            parts = filename.split('_')
            try:
                grid = int([part for part in parts if part.startswith('grid')][0].replace('grid', ''))
                minwin = int([part for part in parts if part.startswith('minwin')][0].replace('minwin', ''))
                maxwin = int([part for part in parts if part.startswith('maxwin')][0].replace('maxwin', ''))
                chr_part = [part for part in parts if part.startswith('chr')][0]
            except (ValueError, IndexError) as e:
                print(f"Skipping file {file} due to error: {e}")
                continue
            
            key = (minwin, maxwin)  
            
            # Read the corresponding info file to get the number of sites
            info_file = file.replace('Report', 'Info')
            number_of_sites = int(next((line.split("Sites:")[1].strip() for line in open(info_file) if "Sites:" in line), None))
            
            results = pd.read_csv(file, comment="/", sep="\t", header=None)
            results.columns = ['Position', 'Likelihood']
            results['grid'] = grid
            results['minwin'] = minwin
            results['maxwin'] = maxwin
            results['model'] = model
            results['sites'] = number_of_sites
            results['chromosome'] = chr_part  
            
            if key not in model_dfs:
                model_dfs[key] = []
            
            model_dfs[key].append(results)
        
        # Concatenate dataframes
        for key, df_list in model_dfs.items():
            combined_df = pd.concat(df_list, ignore_index=True)
            df_dict[f"{model}_minwin{key[0]}_maxwin{key[1]}"] = combined_df
    
    return df_dict

directory = "./"
results = read_OmegaPlus_files(directory)

# Now results is a dictionary where each key corresponds to a specific combination of model, minwin, and maxwin
for key, df in results.items():
    print(f"Combined DataFrame for {key}:\n", df.head())
    # Save to CSV
    df.to_csv(f"{directory}/combined_{key}_grid200kb.csv", index=False)
