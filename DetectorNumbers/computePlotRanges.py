import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

# Read gas densities
gas_data = pd.read_csv('GasDensities.csv')
gas_data.columns = gas_data.columns.str.strip()  # Strip whitespace from column names

# Dictionary to map gas names to CSV file names
gas_to_file = {
    'CF4': 'CF4_Stopping_power_electrons.csv',
    'HeCF4(60% He 40% CF4)': 'HeCF4_Stopping_power_electrons.csv',
    'HeCF4CH4(58% He 37% CF4 5% CH4)': 'HeCF4CH4_Stopping_power_electrons.csv'
}

fig, ax = plt.subplots(figsize=(10, 6))

# Collect energy thresholds table
summary_rows = []


def energy_at_range(target_mm, range_mm, energy_mev):
    """Interpolate energy (MeV) at a given range (mm). Returns NaN if out of bounds."""
    order = np.argsort(range_mm)
    r_sorted = range_mm[order]
    e_sorted = energy_mev[order]

    # Remove duplicate range points for stable interpolation
    r_unique, idx = np.unique(r_sorted, return_index=True)
    e_unique = e_sorted[idx]

    if target_mm < r_unique[0] or target_mm > r_unique[-1]:
        return np.nan

    return np.interp(target_mm, r_unique, e_unique)


for idx, row in gas_data.iterrows():
    gas_name = row['Gas'].strip()
    density = row['density @ 293 K (g/cm3)']
    
    # Get the corresponding stopping power file
    if gas_name not in gas_to_file:
        print(f"Warning: No stopping power file found for {gas_name}")
        continue
    
    filename = gas_to_file[gas_name]
    
    # Read stopping power data, skipping header lines
    try:
        sp_data = pd.read_csv(filename, skiprows=4)
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        continue
    
    # Extract kinetic energy and collision stopping power
    energy_mev = sp_data.iloc[:, 0].values  # Kinetic Energy in MeV
    collision_sp = sp_data.iloc[:, 1].values  # Collision Stp. Pow. in MeV cm2/g
    
    # Extrapolate stopping power down to zero energy
    # Assume 1/E behavior: S(E) ~ 1/E at low energy
    e_min = energy_mev[0]
    sp_min = collision_sp[0]
    # Add low-energy point: extrapolate with 1/E behavior
    extrapolated_energy = np.concatenate([[e_min/10], energy_mev])
    extrapolated_sp = np.concatenate([[sp_min * 10], collision_sp])  # scales as 1/E
    
    energy_mev = extrapolated_energy
    collision_sp = extrapolated_sp
    
    # Calculate linear stopping power: dE/dx = collision_sp * density
    linear_sp = collision_sp * density  # MeV/cm
    
    # Calculate range using constant slowing down approximation
    # R = integral of [1 / (dE/dx)] dE = integral of [1 / (collision_sp * density)] dE
    inverse_sp = 1.0 / linear_sp  # cm/MeV
    
    # Integrate to get range (cumulative from 0 to E)
    range_cm = cumtrapz(inverse_sp, energy_mev, initial=0)
    range_mm = range_cm * 10  # Convert cm to mm

    # --- Added: extract energy at 1 mm and 1 m ranges ---
    e_at_1mm_mev = energy_at_range(1.0, range_mm, energy_mev)
    e_at_1m_mev = energy_at_range(1000.0, range_mm, energy_mev)

    summary_rows.append({
        'Gas': gas_name,
        'density_g_cm3': density,
        'energy_min_at_1mm_keV': e_at_1mm_mev * 1e3 if np.isfinite(e_at_1mm_mev) else np.nan,
        'energy_max_at_1m_keV': e_at_1m_mev * 1e3 if np.isfinite(e_at_1m_mev) else np.nan
    })
    
    # Plot
    label = f"{gas_name} ({density} g/cm³)"
    ax.plot(energy_mev * 1e3, range_mm, label=label)

ax.set_xlabel('Kinetic Energy (keV)')
ax.set_ylabel('Range (mm)')
ax.set_title('Electron Range vs Energy')
ax.set_xlim(left=1, right=10000)  # 1 keV to 10 MeV in keV
ax.axhline(y=1, color='red', linestyle='--', alpha=0.5, linewidth=1)
ax.text(1.5, 1.1, '1 mm', fontsize=10, color='red')
ax.axhline(y=3, color='green', linestyle='--', alpha=0.5, linewidth=1)
ax.text(1.5, 3.1, '3 mm', fontsize=10, color='green')
ax.axhline(y=1000, color='blue', linestyle='--', alpha=0.5, linewidth=1)
ax.text(1.5, 1000.5, '1 m', fontsize=10, color='blue')
ax.set_ylim(bottom=0.01, top=10000)  # 10 microns to 10 meter in mm
ax.legend(title='Gas Mixture', fontsize=8, title_fontsize=9, loc='lower right')
ax.grid(True, alpha=0.3)
ax.set_xscale('log')
ax.set_yscale('log')

plt.tight_layout()
plt.savefig('electron_ranges.png', dpi=150)
# plt.show()

# --- Added: output summary table ---
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv('electron_range_energy_table.csv', index=False)

print("\nEnergy thresholds by gas:")
print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
