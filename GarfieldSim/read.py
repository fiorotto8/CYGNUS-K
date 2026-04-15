import argparse
import csv
import ctypes
from pathlib import Path

import numpy as np
import ROOT
import Garfield


def electron_diffusion_at(gas, e_field_v_cm):
	dl = ctypes.c_double(0.0)
	dt = ctypes.c_double(0.0)
	ok = gas.ElectronDiffusion(e_field_v_cm, 0.0, 0.0, 0.0, 0.0, 0.0, dl, dt)
	return ok, dl.value, dt.value


def electron_transport_at(gas, e_field_v_cm):
	ok_diff, dl, dt = electron_diffusion_at(gas, e_field_v_cm)

	alpha = ctypes.c_double(0.0)
	ok_townsend = gas.ElectronTownsend(
		e_field_v_cm, 0.0, 0.0, 0.0, 0.0, 0.0, alpha
	)

	eta = ctypes.c_double(0.0)
	ok_attachment = gas.ElectronAttachment(
		e_field_v_cm, 0.0, 0.0, 0.0, 0.0, 0.0, eta
	)

	vx = ctypes.c_double(0.0)
	vy = ctypes.c_double(0.0)
	vz = ctypes.c_double(0.0)
	ok_velocity = gas.ElectronVelocity(
		e_field_v_cm, 0.0, 0.0, 0.0, 0.0, 0.0, vx, vy, vz
	)

	speed = (vx.value ** 2 + vy.value ** 2 + vz.value ** 2) ** 0.5 if ok_velocity else float("nan")

	return {
		"ok_diff": ok_diff,
		"ok_townsend": ok_townsend,
		"ok_attachment": ok_attachment,
		"ok_velocity": ok_velocity,
		"dl": dl,
		"dt": dt,
		"townsend": alpha.value,
		"attachment": eta.value,
		"vx": vx.value,
		"vy": vy.value,
		"vz": vz.value,
		"speed": speed,
	}


def refine_minimum_with_quadratic(e_vals, y_vals):
	"""Refine minimum using a 3-point quadratic fit in log10(E)."""
	n = len(e_vals)
	if n == 0:
		return None, None, False
	if n < 3:
		i = int(np.argmin(y_vals))
		return e_vals[i], y_vals[i], False

	i0 = int(np.argmin(y_vals))
	if i0 == 0 or i0 == n - 1:
		return e_vals[i0], y_vals[i0], False

	x = np.log10(np.array(e_vals[i0 - 1:i0 + 2], dtype=float))
	y = np.array(y_vals[i0 - 1:i0 + 2], dtype=float)

	a, b, c = np.polyfit(x, y, 2)
	if abs(a) < 1e-18 or a <= 0:
		return e_vals[i0], y_vals[i0], False

	xv = -b / (2.0 * a)
	if xv < x.min() or xv > x.max():
		return e_vals[i0], y_vals[i0], False

	e_refined = float(10.0 ** xv)
	y_refined = float(a * xv * xv + b * xv + c)
	return e_refined, y_refined, True


def main():
	parser = argparse.ArgumentParser(
		description="Read Garfield gas table, plot transport quantities, and report diffusion values."
	)
	parser.add_argument(
		"gasfile",
		help="Path to the .gas file (for example: tables/heCf4_60-40_1000mbar20C.gas)",
	)
	args = parser.parse_args()

	gas_file = Path(args.gasfile)
	if not gas_file.exists():
		raise FileNotFoundError(f"Gas file not found: {gas_file}")

	table_name = gas_file.stem
	out_dir = Path("plots")
	out_dir.mkdir(parents=True, exist_ok=True)

	gas = ROOT.Garfield.MediumMagboltz()
	if not gas.LoadGasFile(str(gas_file)):
		raise RuntimeError(f"Could not load gas file: {gas_file}")

	gas.PrintGas()

	view = ROOT.Garfield.ViewMedium(gas)

	c_v = ROOT.TCanvas("cV", "", 700, 600)
	view.SetCanvas(c_v)
	view.PlotElectronVelocity()
	c_v.SaveAs(str(out_dir / f"{table_name}_electron_velocity.png"))

	c_d = ROOT.TCanvas("cD", "", 700, 600)
	view.SetCanvas(c_d)
	view.PlotElectronDiffusion()
	c_d.SaveAs(str(out_dir / f"{table_name}_electron_diffusion.png"))

	c_t = ROOT.TCanvas("cT", "", 700, 600)
	view.SetCanvas(c_t)
	view.PlotElectronTownsend()
	c_t.SaveAs(str(out_dir / f"{table_name}_electron_townsend.png"))

	c_a = ROOT.TCanvas("cA", "", 700, 600)
	view.SetCanvas(c_a)
	view.PlotElectronAttachment()
	c_a.SaveAs(str(out_dir / f"{table_name}_electron_attachment.png"))

	efields = ROOT.std.vector["double"]()
	bfields = ROOT.std.vector["double"]()
	angles = ROOT.std.vector["double"]()
	gas.GetFieldGrid(efields, bfields, angles)

	if efields.size() == 0:
		print("No electric-field grid found in gas table.")
		return

	# Diffusion at fixed reference fields.
	e_query_1kv = 1000.0
	e_query = 2000.0
	ok_query_1kv, dl_query_1kv, dt_query_1kv = electron_diffusion_at(gas, e_query_1kv)
	if ok_query_1kv:
		print(f"At E = {e_query_1kv:.1f} V/cm (1 kV/cm), B = 0 T:")
		print(f"  Longitudinal diffusion (DL) = {dl_query_1kv:.6g}")
		print(f"  Transverse diffusion (DT)   = {dt_query_1kv:.6g}")
	else:
		print(f"Could not evaluate electron diffusion at E = {e_query_1kv:.1f} V/cm.")

	ok_query, dl_query, dt_query = electron_diffusion_at(gas, e_query)
	if ok_query:
		print(f"At E = {e_query:.1f} V/cm (2 kV/cm), B = 0 T:")
		print(f"  Longitudinal diffusion (DL) = {dl_query:.6g}")
		print(f"  Transverse diffusion (DT)   = {dt_query:.6g}")
	else:
		print(f"Could not evaluate electron diffusion at E = {e_query:.1f} V/cm.")

	# Find minima on the tabulated electric-field grid and refine by interpolation.
	e_grid = []
	dl_grid = []
	dt_grid = []

	for i in range(efields.size()):
		e_val = efields[i]
		ok, dl_val, dt_val = electron_diffusion_at(gas, e_val)
		if not ok:
			continue
		e_grid.append(float(e_val))
		dl_grid.append(float(dl_val))
		dt_grid.append(float(dt_val))

	e_at_min_dl, min_dl, dl_interpolated = refine_minimum_with_quadratic(e_grid, dl_grid)
	e_at_min_dt, min_dt, dt_interpolated = refine_minimum_with_quadratic(e_grid, dt_grid)

	if e_at_min_dl is not None:
		method = "interpolated" if dl_interpolated else "grid"
		print(f"Minimum longitudinal diffusion ({method}):")
		print(f"  DL_min = {min_dl:.6g} at E = {e_at_min_dl:.6g} V/cm")
	else:
		print("Could not determine minimum longitudinal diffusion on grid.")

	if e_at_min_dt is not None:
		method = "interpolated" if dt_interpolated else "grid"
		print(f"Minimum transverse diffusion ({method}):")
		print(f"  DT_min = {min_dt:.6g} at E = {e_at_min_dt:.6g} V/cm")
	else:
		print("Could not determine minimum transverse diffusion on grid.")

	if e_at_min_dl is not None:
		rows = []
		for label, e_val in [
			("E_1kV_per_cm", e_query_1kv),
			("E_2kV_per_cm", e_query),
			("E_at_DL_min", e_at_min_dl),
		]:
			vals = electron_transport_at(gas, e_val)
			dl_out = vals["dl"] if label != "E_at_DL_min" else min_dl
			rows.append({
				"point": label,
				"electric_field(V/cm)": e_val,
				"DL(sqrt(cm))": dl_out,
				"DT(sqrt(cm))": vals["dt"],
				"Townsend": vals["townsend"],
				"Attachment(cm-1)": vals["attachment"],
				"Velocity(cm/ns)": vals["speed"],
				"vx": vals["vx"],
				"vy": vals["vy"],
				"vz": vals["vz"],
				"ok_diff": vals["ok_diff"],
				"ok_townsend": vals["ok_townsend"],
				"ok_attachment": vals["ok_attachment"],
				"ok_velocity": vals["ok_velocity"],
			})

		csv_path = out_dir / f"{table_name}_transport_points.csv"
		with csv_path.open("w", newline="") as f:
			writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
			writer.writeheader()
			writer.writerows(rows)
		print(f"Saved transport values to {csv_path}")


if __name__ == "__main__":
	main()