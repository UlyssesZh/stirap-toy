import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_qutip_csv(file_name):
    # Standard loading function for the provided QuTiP CSV format
    df = pd.read_csv(file_name, comment='#', header=None)
    def complex_to_real(val):
        if isinstance(val, str):
            # Remove spaces and convert complex string to float (real part)
            return float(complex(val.replace(' ', '')).real)
        return val
    for col in df.columns:
        df[col] = df[col].apply(complex_to_real)
    df.columns = ['time', 'n1', 'n2', 'pump', 'stokes']
    return df

# Load the new datasets
fstirap_new = load_qutip_csv('fstirap_0+1.csv')
rfstirap_new = load_qutip_csv('rfstirap_0+1.csv')

# 1. Mirror the pulses for fstirap section
fstirap_mirrored = fstirap_new.copy()
fstirap_mirrored['pump'] = fstirap_new['pump'].values[::-1]
fstirap_mirrored['stokes'] = fstirap_new['stokes'].values[::-1]

# 2. Shift rfstirap time by 4.0 ms to start at 2.0 ms (original -2.0 + 4.0 = 2.0)
rfstirap_shifted = rfstirap_new.copy()
rfstirap_shifted['time'] = rfstirap_shifted['time'] + 4.0

# 3. Merge datasets
merged_df_new = pd.concat([fstirap_mirrored, rfstirap_shifted], ignore_index=True)

# 4. Save to CSV
merged_df_new.to_csv('merged_stirap_0+1.csv', index=False)

import matplotlib.pyplot as plt

# 5. Plotting
fig, ax1 = plt.subplots(figsize=(10, 6))

# --- Primary Y-Axis (Populations) ---
line1, = ax1.plot(merged_df_new['time'], merged_df_new['n1'],
                  label='$n_1$', color='red', linewidth=2)
line2, = ax1.plot(merged_df_new['time'], merged_df_new['n2'],
                  label='$n_2$', color='blue', linewidth=2)

ax1.set_xlabel('Time (ms)', fontsize=20)
ax1.set_ylabel('Probability Amplitude', fontsize=20)
ax1.tick_params(axis='y', labelsize=18)
ax1.tick_params(axis='x', labelsize=18)
ax1.set_xlim(-2.0, 6.0)
ax1.set_ylim(bottom=0.0, top=1.01)
# --- Secondary Y-Axis (Pulses) ---
ax2 = ax1.twinx()  # Create the second axis
line3, = ax2.plot(merged_df_new['time'], 2000*merged_df_new['pump'],
                  label='Pump', linestyle=':', color='red', alpha=0.7)
line4, = ax2.plot(merged_df_new['time'], 2000*merged_df_new['stokes'],
                  label='Stokes', linestyle=':', color='blue', alpha=0.7)

ax2.set_ylabel('Pulse Amplitude', fontsize=20)
ax2.tick_params(axis='y', labelsize=18)
ax2.set_ylim(bottom=0.0, top=4000)
# --- Combined Legend ---
# Since we have two axes, we gather labels from both to make one legend
lines = [line1, line2, line3, line4]
labels = [l.get_label() for l in lines]


plt.tight_layout()

# Save image
plt.savefig('merged_stirap_0+1_final_plot.pdf')
plt.show()
