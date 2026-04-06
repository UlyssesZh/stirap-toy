import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_qutip_csv(file_name):
    df = pd.read_csv(file_name, comment='#', header=None)
    def complex_to_real(val):
        if isinstance(val, str):
            return float(complex(val.replace(' ', '')).real)
        return val
    for col in df.columns:
        df[col] = df[col].apply(complex_to_real)
    df.columns = ['time', 'n1', 'n2', 'pump', 'stokes']
    return df

fstirap = load_qutip_csv('fstirap_10mK.csv')
rfstirap = load_qutip_csv('rfstirap_10mK.csv')

# Mirroring the pulses for fstirap
fstirap_mirrored = fstirap.copy()
fstirap_mirrored['pump'] = fstirap['pump'].values[::-1]
fstirap_mirrored['stokes'] = fstirap['stokes'].values[::-1]

# Shift rfstirap time by 4.0 to start at 2.0
rfstirap_shifted = rfstirap.copy()
rfstirap_shifted['time'] = rfstirap_shifted['time'] + 4.0

# Merge datasets
merged_df = pd.concat([fstirap_mirrored, rfstirap_shifted], ignore_index=True)

# Plotting
plt.figure(figsize=(10, 6))

# Plot Populations
plt.plot(merged_df['time'], merged_df['n1'], label='$n_1$', color='red', linewidth=2)
plt.plot(merged_df['time'], merged_df['n2'], label='$n_2$', color='blue', linewidth=2)

# Plot Pulses
plt.plot(merged_df['time'], merged_df['pump'], label='Pump', linestyle=':', color='red', alpha=0.7)
plt.plot(merged_df['time'], merged_df['stokes'], label='Stokes', linestyle=':', color='blue', alpha=0.7)

# Vertical axis title
plt.ylabel(r'$\langle n \rangle$', fontsize=14)
plt.xlabel('Time (ms)', fontsize=12)

# Grid
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Place legend at the middle top
# Using ncol=4 to arrange labels horizontally and avoid overlapping with data
plt.legend(loc='upper center', ncol=4, frameon=True)

# Adjust layout to make sure legend fits if needed,
# but loc='upper center' is inside the axis.
plt.tight_layout()

# Save image
plt.savefig('stirap_legend_top_center.png')
