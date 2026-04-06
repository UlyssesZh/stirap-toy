import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- 1. FUNCTION DEFINITIONS ---
def load_qutip_csv(file_name):
    try:
        df = pd.read_csv(file_name, comment='#', header=None)
        def complex_to_real(val):
            if isinstance(val, str):
                return float(complex(val.replace(' ', '')).real)
            return val
        for col in df.columns:
            df[col] = df[col].apply(complex_to_real)
        df.columns = ['time', 'n1', 'n2', 'pump', 'stokes']
        return df
    except FileNotFoundError:
        print(f"File {file_name} not found.")
        return None

def scale_to_custom_range(df, start, end):
    """Normalizes time and scales it to a specific [start, end] range."""
    t = df['time'].values
    t_norm = (t - t.min()) / (t.max() - t.min())
    df['time'] = t_norm * (end - start) + start
    return df

# --- 2. DATA LOADING ---
f_raw = load_qutip_csv('fstirap_0+1_10mk.csv')
rf_raw = load_qutip_csv('rfstirap_0+1_10mk.csv')

try:
    data_f = np.load('negativity_evolution_fstirap.npz')
    data_rf = np.load('negativity_evolution_rfstirap.npz')
    neg_f = data_f['negativity']
    neg_rf = data_rf['negativity']
except Exception as e:
    print(f"Error loading negativity files: {e}")
    neg_f, neg_rf = np.zeros(len(f_raw)), np.zeros(len(rf_raw))

# --- 3. PROCESSING & SCALING ---

# SEGMENT 1: Let's assume this is a 'pre-buffer' mirrored segment.
# Scaling it to be the first 2ms of the plot (-4.0 to -2.0) or similar.
# For a continuous plot starting at -2.0, we will scale Segment 1 to [-2, 2]
# and Segment 2 to [2, 6] to fill the 8ms window.

f_mirrored = scale_to_custom_range(f_raw.copy(), -2.0, 2.0)
f_mirrored['pump'] = f_mirrored['pump'].values[::-1]
f_mirrored['stokes'] = f_mirrored['stokes'].values[::-1]
neg_f_processed = neg_f

# SEGMENT 2: The standard sequence filling the second half [2.0 to 6.0]
rf_seg = scale_to_custom_range(rf_raw.copy(), 2.0, 6.0)
neg_rf_processed = neg_rf

# --- 4. MERGING ---
merged_df = pd.concat([f_mirrored, rf_seg], ignore_index=True)
merged_negativity = np.concatenate([neg_f_processed, neg_rf_processed])

# --- 5. PLOTTING ---
fig, ax1 = plt.subplots(figsize=(10, 6))

# Primary Axis
ax1.plot(merged_df['time'], merged_df['n1'], label='$n_1$', color='red', lw=2)
ax1.plot(merged_df['time'], merged_df['n2'], label='$n_2$', color='blue', lw=2)
ax1.plot(merged_df['time'], merged_negativity, label=r'Negativity $\mathcal{N}$',
                 color='black', ls='--', lw=2.5)

ax1.set_xlabel('Time (ms)', fontsize=20)
ax1.set_ylabel('Probability / Negativity', fontsize=20)
ax1.tick_params(axis='both', labelsize=20)
ax1.set_xlim(-2.0, 6.0) # Set the requested window
ax1.set_ylim(0, 1.01)

# Secondary Axis: Pulses
ax2 = ax1.twinx()
ax2.plot(merged_df['time'], 2000*merged_df['pump'], label='Pump Pulse', ls=':', color='red', lw=1.5)
ax2.plot(merged_df['time'], 2000*merged_df['stokes'], label='Stokes Pulse', ls=':', color='blue', lw=1.5)
ax2.set_ylabel('Pulse Amplitude', fontsize=20)
ax2.set_ylim(0, 4000)
ax2.tick_params(axis='y', labelsize=20)

# Combined Legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

plt.tight_layout()
plt.savefig('merged_fstirap_8ms_plot.pdf')
plt.show()
