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

def scale_to_1ms(df):
    """Normalizes time to 0.0 - 1.0 ms range."""
    t = df['time'].values
    df['time'] = (t - t.min()) / (t.max() - t.min())
    return df

# --- 2. DATA LOADING ---
f_raw = load_qutip_csv('fstirap_0+1.csv')
rf_raw = load_qutip_csv('rfstirap_0+1.csv')

# Load Coherence data (modulus from the density matrix)
try:
    data_f = np.load('coherence_evolution_fstirap.npz')
    data_rf = np.load('coherence_evolution_rfstirap.npz')
    coh_f = data_f['coherence']
    coh_rf = data_rf['coherence']
except Exception as e:
    print(f"Error loading coherence files: {e}")
    coh_f, coh_rf = np.zeros(len(f_raw)), np.zeros(len(rf_raw))

# --- 3. PROCESSING & PULSE MIRRORING ---

# SEGMENT 1 (0 to 1 ms): Mirror the pulses ONLY
f_mirrored = scale_to_1ms(f_raw.copy())
# Reverse only the pulse columns to create a "reversed" sequence
f_mirrored['pump'] = f_mirrored['pump'].values[::-1]
f_mirrored['stokes'] = f_mirrored['stokes'].values[::-1]
# Note: Populations (n1, n2) are NOT reversed here to keep the
# continuity with Segment 2, assuming Segment 1 ends where 2 begins.
coh_f_processed = coh_f # Coherence remains chronological

# SEGMENT 2 (1 to 2 ms): Standard sequence
rf_seg = scale_to_1ms(rf_raw.copy())
rf_seg['time'] = rf_seg['time'] + 1.0
coh_rf_processed = coh_rf

# --- 4. MERGING ---
merged_df = pd.concat([f_mirrored, rf_seg], ignore_index=True)
merged_coherence = np.concatenate([coh_f_processed, coh_rf_processed])

# --- 5. PLOTTING ---
fig, ax1 = plt.subplots(figsize=(10, 6))

# Primary Axis: Populations and Coherence
line1, = ax1.plot(merged_df['time'], merged_df['n1'], label='$n_1$', color='red', lw=2)
line2, = ax1.plot(merged_df['time'], merged_df['n2'], label='$n_2$', color='blue', lw=2)
line_c, = ax1.plot(merged_df['time'], merged_coherence, label=r'$|\rho_{coh}|$',
                  color='black', ls='--', lw=2.5)

ax1.set_xlabel('Time (ms)', fontsize=20)
ax1.set_ylabel('Probability / Coherence Modulus', fontsize=20)
ax1.tick_params(axis='y', labelsize=18)
ax1.tick_params(axis='x', labelsize=18)
ax1.set_xlim(0, 2.0)
ax1.set_ylim(0, 1.0) # Based on your image scaleax1.grid(True, linestyle=':', alpha=0.6)

# Secondary Axis: Pulses
ax2 = ax1.twinx()
# Pulse scaling 2000 to match your axis
line3, = ax2.plot(merged_df['time'], 2000*merged_df['pump'], label='Pump Pulse',
                  ls=':', color='red', alpha=0.4)
line4, = ax2.plot(merged_df['time'], 2000*merged_df['stokes'], label='Stokes Pulse',
                  ls=':', color='blue', alpha=0.4)
ax2.set_ylabel('Pulse Amplitude', fontsize=20)
ax2.set_ylim(0, 4500)
ax2.tick_params(axis='y', labelsize=18)
plt.tight_layout()
plt.savefig('merged_fstirap_1K_final_plot.pdf')
plt.show()
