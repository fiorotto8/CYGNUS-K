from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DRIFT_LENGTH_CM = 50.0


def find_input_dir(repo_root: Path) -> Path:
	candidates = [
		repo_root / "GarfieldDiffusion" / "plots",
		repo_root / "GarfieldSim" / "plots",
	]
	for directory in candidates:
		if directory.exists() and directory.is_dir():
			return directory
	raise FileNotFoundError(
		"Could not find input folder. Expected one of: GarfieldDiffusion/plots or GarfieldSim/plots"
	)


def compute_two_sigma(diffusion_coeff_sqrt_cm: float, drift_length_cm: float) -> float:
	# With D in sqrt(cm), sigma = D * sqrt(L), so 2sigma = 2 * D * sqrt(L).
	return 2.0 * diffusion_coeff_sqrt_cm * np.sqrt(drift_length_cm)


def make_comparison_plots(df: pd.DataFrame, output_png: Path) -> None:
	fig, (ax_dl, ax_dt) = plt.subplots(1, 2, figsize=(12, 5), sharex=True)

	for gas_name, gdf in df.groupby("gas"):
		gdf_sorted = gdf.sort_values("electric_field(V/cm)")
		x = gdf_sorted["electric_field(V/cm)"].to_numpy()
		ax_dl.plot(x, gdf_sorted["DL_2sigma(mm)"].to_numpy(), marker="o", label=gas_name)
		ax_dt.plot(x, gdf_sorted["DT_2sigma(mm)"].to_numpy(), marker="o", label=gas_name)

	ax_dl.set_title("Longitudinal diffusion (2σ)")
	ax_dt.set_title("Transverse diffusion (2σ)")
	ax_dl.set_xlabel("Electric field (V/cm)")
	ax_dt.set_xlabel("Electric field (V/cm)")
	ax_dl.set_ylabel("2σ diffusion over 50 cm (mm)")
	ax_dt.set_ylabel("2σ diffusion over 50 cm (mm)")
	ax_dl.grid(True, alpha=0.3)
	ax_dt.grid(True, alpha=0.3)
	ax_dl.legend(fontsize=8)
	ax_dt.legend(fontsize=8)

	fig.suptitle("Diffusion comparison by gas")
	fig.tight_layout()
	fig.savefig(output_png, dpi=150)
	plt.close(fig)


def main() -> None:
	this_file = Path(__file__).resolve()
	repo_root = this_file.parent.parent
	input_dir = find_input_dir(repo_root)
	output_csv = this_file.parent / "diffusion_2sigma_50cm_summary.csv"
	output_plot = this_file.parent / "diffusion_2sigma_vs_field_by_gas.png"

	csv_files = sorted(input_dir.glob("*_transport_points.csv"))
	if not csv_files:
		raise FileNotFoundError(f"No *_transport_points.csv files found in {input_dir}")

	rows = []
	for csv_file in csv_files:
		gas_description = csv_file.name.replace("_transport_points.csv", "")
		df = pd.read_csv(csv_file)
		df.columns = [c.strip() for c in df.columns]

		required_cols = ["point", "electric_field(V/cm)", "DL(sqrt(cm))", "DT(sqrt(cm))"]
		missing = [c for c in required_cols if c not in df.columns]
		if missing:
			raise ValueError(f"Missing columns in {csv_file.name}: {missing}")

		for _, r in df.iterrows():
			dl = float(r["DL(sqrt(cm))"])
			dt = float(r["DT(sqrt(cm))"])
			rows.append(
				{
					"gas": gas_description,
					"point": r["point"],
					"electric_field(V/cm)": float(r["electric_field(V/cm)"]),
					"drift_length(cm)": DRIFT_LENGTH_CM,
					"DL(sqrt(cm))": dl,
					"DT(sqrt(cm))": dt,
					"DL_2sigma(cm)": compute_two_sigma(dl, DRIFT_LENGTH_CM),
					"DT_2sigma(cm)": compute_two_sigma(dt, DRIFT_LENGTH_CM),
					"DL_2sigma(mm)": 10.0 * compute_two_sigma(dl, DRIFT_LENGTH_CM),
					"DT_2sigma(mm)": 10.0 * compute_two_sigma(dt, DRIFT_LENGTH_CM),
				}
			)

	out_df = pd.DataFrame(rows)
	out_df.to_csv(output_csv, index=False)
	make_comparison_plots(out_df, output_plot)

	print(f"Input directory: {input_dir}")
	print(f"Read {len(csv_files)} files and {len(out_df)} rows.")
	print(f"Saved summary: {output_csv}")
	print(f"Saved plot: {output_plot}")


if __name__ == "__main__":
	main()
