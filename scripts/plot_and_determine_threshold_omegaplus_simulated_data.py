import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import math

# Load dataframes
df_sagit = pd.read_csv("combined_sagit_grid10_minwin10000_maxwin100000.csv")
df_nemo = pd.read_csv("combined_nemo_grid10_minwin7500_maxwin100000.csv")

# Combine the dataframes
combined_df = pd.concat([df_sagit, df_nemo], ignore_index=True)

# Define window size for the additional plot
window_size = 500

# Function to assign window based on position
def assign_window(position):
    return int(position // window_size)

combined_df['Window'] = combined_df['Position'].apply(assign_window)

# Calculate the mean likelihood and standard error for each model in each window
grouped = combined_df.groupby(['model', 'Window']).agg(
    mean_likelihood=('Likelihood', 'mean'),
    std_likelihood=('Likelihood', 'std'),
    count=('Likelihood', 'size')
).reset_index()

# Calculate standard error and confidence intervals
grouped['se'] = grouped['std_likelihood'] / np.sqrt(grouped['count'])
z = norm.ppf(0.995)
grouped['lower'] = grouped['mean_likelihood'] - z * grouped['se']
grouped['upper'] = grouped['mean_likelihood'] + z * grouped['se']

# Calculate 99% quantile of likelihood values for each model
quantiles = combined_df.groupby('model')['Likelihood'].quantile(0.99).reset_index()
quantiles.columns = ['model', '99%_quantile']
print(quantiles)

# Merge quantiles with combined_df to check which points are above the 99% quantile
combined_df = combined_df.merge(quantiles, on='model')
combined_df['above_99_quantile'] = combined_df['Likelihood'] > combined_df['99%_quantile']
num_models = combined_df['model'].nunique()
models = combined_df['model'].unique()

# Create subplots in a 1x2 layout
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=False)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Loop through each model and plot in its own subplot
for i, (ax, model) in enumerate(zip(axes, models)):
    model_data = grouped[grouped['model'] == model]
    ax.plot(model_data['Window'] * window_size, model_data['mean_likelihood'], label=model, color=colors[i % len(colors)])
    ax.fill_between(model_data['Window'] * window_size, model_data['lower'], model_data['upper'], color=colors[i % len(colors)], alpha=0.3)
    
    # Add horizontal line for 99% quantile
    quantile_value = quantiles[quantiles['model'] == model]['99%_quantile'].values[0]
    ax.axhline(y=quantile_value, color='red', linestyle='--', label='99% Quantile')
    ax.set_ylim(0, 25)
    ax.set_title(f'Likelihood for {model}')
    ax.set_xlabel('Position')
    ax.set_ylabel('Mean Likelihood')
    # Add quantile information as text annotation in the top right corner
    ax.text(0.95, 0.95, f'99% Quantile: {quantile_value:.2f}', transform=ax.transAxes, 
            verticalalignment='top', horizontalalignment='right', 
            bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))
    ax.legend(loc='lower left')
plt.tight_layout()
output_plot_path = "likelihood_plot_grid10.png"
plt.savefig(output_plot_path)
print(f"Plot saved to {output_plot_path}")
plt.show()

# Create subplots for the additional scatter plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=False)

# Loop through each model and plot in its own subplot with dots
for i, (ax, model) in enumerate(zip(axes, models)):
    model_data = grouped[grouped['model'] == model]
    ax.scatter(model_data['Window'] * window_size, model_data['mean_likelihood'], s=4, label=model, color=colors[i % len(colors)])
    ax.fill_between(model_data['Window'] * window_size, model_data['lower'], model_data['upper'], color=colors[i % len(colors)], alpha=0.3)
    
    # Add horizontal line for 99% quantile
    quantile_value = quantiles[quantiles['model'] == model]['99%_quantile'].values[0]
    ax.axhline(y=quantile_value, color='red', linestyle='--', label='99% Quantile')
    ax.set_ylim(0, 25)
    ax.set_title(f'Likelihood for {model}')
    ax.set_xlabel('Position')
    ax.set_ylabel('Mean Likelihood')
    
    # Add quantile information as text annotation in the top right corner
    ax.text(0.95, 0.95, f'99% Quantile: {quantile_value:.2f}', transform=ax.transAxes, 
            verticalalignment='top', horizontalalignment='right', 
            bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))
    ax.legend(loc='lower left')
plt.tight_layout()

output_scatter_plot_path = "likelihood_scatter_plot_grid10.png"
plt.savefig(output_scatter_plot_path)
print(f"Scatter plot saved to {output_scatter_plot_path}")
plt.show()

print("Data points above 99% quantile:")
print(combined_df[combined_df['above_99_quantile']])
