import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

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

def scale_to_1ms(df):
    """Normalizes time to 0.0 - 1.0 ms range."""
    t = df['time'].values
    df['time'] = (t - t.min()) / (t.max() - t.min())
    return df

# --- 2. DATA LOADING ---
f_raw = load_qutip_csv('fstirap_0+1.csv')
rf_raw = load_qutip_csv('rfstirap_0+1.csv')

# Load negativity data (modulus from the density matrix)
try:
    data_f = np.load('negativity_evolution_fstirap.npz')
    data_rf = np.load('negativity_evolution_rfstirap.npz')
    neg_f = data_f['negativity']
    neg_rf = data_rf['negativity']
except Exception as e:
    print(f"Error loading negativity files: {e}")
    neg_f, neg_rf = np.zeros(len(f_raw)), np.zeros(len(rf_raw))

# --- 3. PROCESSING & PULSE MIRRORING ---

# SEGMENT 1 (0 to 1 ms): Mirror the pulses ONLY
f_mirrored = scale_to_1ms(f_raw.copy())
# Reverse only the pulse columns to create a "reversed" sequence
f_mirrored['pump'] = f_mirrored['pump'].values[::-1]
f_mirrored['stokes'] = f_mirrored['stokes'].values[::-1]
# Note: Populations (n1, n2) are NOT reversed here to keep the
# continuity with Segment 2, assuming Segment 1 ends where 2 begins.
neg_f_processed = neg_f # negativity remains chronological

# SEGMENT 2 (1 to 2 ms): Standard sequence
rf_seg = scale_to_1ms(rf_raw.copy())
rf_seg['time'] = rf_seg['time'] + 1.0
neg_rf_processed = neg_rf

# --- 4. MERGING ---
merged_df = pd.concat([f_mirrored, rf_seg], ignore_index=True)
merged_negativity = np.concatenate([neg_f_processed, neg_rf_processed])

# --- 5. PLOTTING ---
fig, ax1 = plt.subplots(figsize=(10, 6))

# Primary Axis: Populations and negativity
line1, = ax1.plot(merged_df['time'], merged_df['n1'], label='$n_1$', color='red', lw=2)
line2, = ax1.plot(merged_df['time'], merged_df['n2'], label='$n_2$', color='blue', lw=2)
line_c, = ax1.plot(merged_df['time'], merged_negativity, label=r'$|\rho_{neg}|$',
                  color='black', ls='--', lw=2.5)

# --- KEY CHANGES FOR ZERO OVERLAP ---


ticks = ax1.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax1.set_xticklabels(['', '0.5', '1.0', '1.5', ''])
# ------------------------------------

# 2. Style Labels: Match Image 1
ax1.set_xlabel('Time (ms)', fontsize=20, labelpad=-2)
ax1.set_ylabel('Probability / Negativity', fontsize=20)

# 3. Tick Parameters: Inward ticks like the publication style
ax1.tick_params(axis='both', which='major', labelsize=20, direction='in', length=6)
ax1.set_xlim(0, 2.0)
ax1.set_ylim(0, 0.9)

# Secondary Axis: Pulses
ax2 = ax1.twinx()
line3, = ax2.plot(merged_df['time'], 2000*merged_df['pump'], ls=':', color='red', lw=1.5)
line4, = ax2.plot(merged_df['time'], 2000*merged_df['stokes'], ls=':', color='blue', lw=1.5)

ax2.set_ylabel('Pulse Amplitude', fontsize=20)
ax2.set_ylim(0, 4000)

# Apply inward ticks to secondary axis
ax2.tick_params(axis='y', labelsize=20, direction='in', length=6)

plt.tight_layout()
plt.savefig('merged_fstirap_1K_final_plot.pdf')
plt.show()
