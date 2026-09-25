"""Bounded diagnostics for the TOI-1231 b gas-giant reproducer.

Run from this directory with: conda run -n photochem python -u diagnose_initial.py
"""

import argparse
import json
import os
import tempfile
import time

import numpy as np
import yaml

import test as case
from photochem.extensions import gasgiants as current


def profile_error(model):
    pressure = model.wrk.pressure_hydro
    desired_pressure = model.gdat.P_desired[::-1]
    desired_temperature = model.gdat.T_desired[::-1]
    desired_edd = model.gdat.Kzz_desired[::-1]
    logp = np.log10(pressure)
    expected_temperature = np.interp(logp, np.log10(desired_pressure), desired_temperature)
    expected_edd = np.interp(logp, np.log10(desired_pressure), np.log10(desired_edd))
    return (
        float(np.max(np.abs(model.var.temperature - expected_temperature))),
        float(np.max(np.abs(np.log10(model.var.edd) - expected_edd))),
    )


def snapshot(model, label, step=None, converged=False, give_up=False):
    pressure = model.wrk.pressure_hydro
    temperature_error, edd_error = profile_error(model)
    history = model.wrk.t_history
    started = label != "initialized"
    out = {
        "label": label,
        "step": step,
        "time": float(model.wrk.tn),
        "dt_accepted": float(history[0] - history[1]) if step else None,
        "nsteps_segment": int(model.wrk.nsteps),
        "nsteps_total": int(model.wrk.nsteps_total) if started else None,
        "longdy": float(model.wrk.longdy),
        "longdydt": float(model.wrk.longdydt),
        "toa_pressure": float(pressure[-1]),
        "bottom_pressure": float(pressure[0]),
        "max_temperature_error": temperature_error,
        "max_log10_edd_error": edd_error,
        "toa_updates": int(model.wrk.n_toa_pressure_updates),
        "errors_total": int(model.wrk.nerrors_total) if started else None,
        "converged": converged,
        "give_up": give_up,
    }
    print(json.dumps(out, allow_nan=True), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", choices=("old", "new"), default="new")
    parser.add_argument("--nonpersistent", action="store_true")
    parser.add_argument("--clear-persistent", action="store_true")
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--interval", type=int, default=50)
    parser.add_argument("--reinit-at", type=int)
    parser.add_argument("--detail-start", type=int, default=0)
    parser.add_argument("--detail-end", type=int, default=0)
    parser.add_argument("--verbose-start", type=int, default=0)
    parser.add_argument("--verbose-end", type=int, default=0)
    parser.add_argument("--verbose-level", type=int, default=1)
    parser.add_argument("--probe-at", type=int)
    parser.add_argument("--continuous-c3h6", action="store_true")
    args = parser.parse_args()
    temporary_mechanism = None
    old_module = case.gasgiants_0_6_7
    if args.version == "new":
        if args.continuous_c3h6:
            with open("photochem_rxns.yaml", encoding="utf-8") as stream:
                mechanism_text = stream.read()
            original = "[173.78, 13.606, -2.1834, 0.1228, -45.669, -27.1, 351.0]"
            adjusted = "[173.78, 13.606, -2.1834, 0.1228, -45.669, -133.53375344651629, 351.0]"
            if mechanism_text.count(original) != 1:
                raise RuntimeError("Expected one C3H6 high-temperature coefficient row")
            with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as stream:
                stream.write(mechanism_text.replace(original, adjusted))
                temporary_mechanism = stream.name

        class DiagnosticGasGiant(current.EvoAtmosphereGasGiant):
            def __init__(self, mechanism_file, *profile_args, **kwargs):
                super().__init__(temporary_mechanism or mechanism_file, *profile_args, **kwargs)

            def initialize_atmosphere_p(self, *profile_args, **kwargs):
                if args.nonpersistent:
                    kwargs["persistent"] = False
                    kwargs["maintain_toa_pressure"] = None
                    kwargs["target_pressure"] = None
                return super().initialize_atmosphere_p(*profile_args, **kwargs)

        class DiagnosticModule:
            EvoAtmosphereGasGiant = DiagnosticGasGiant

        selected_module = DiagnosticModule
    else:
        selected_module = old_module
    # test.py is a reproducer that the user may switch between the two
    # extensions. Force both names so its local selection cannot bypass a
    # diagnostic variant.
    case.gasgiants_0_6_7 = selected_module
    case.gasgiants_0_9_0 = selected_module
    start = time.monotonic()
    model = case.initialize_photochem(
        spectrum="TOI1231b_spectrum.txt",
        planet_mass=case.planets.TOI1231b.mass,
        planet_radius=case.planets.TOI1231b.radius,
        metallicity=80,
        x_acc=None,
        CtoO=1.0,
        N_depletion=None,
        climate_filename="TOI1231b_MH=2.000_CO=1.000_Tint=50.0.pkl",
        Kzz=1e7,
        initial_cond_with_quenching=True,
    )
    if args.version == "new" and not isinstance(model, DiagnosticGasGiant):
        raise RuntimeError("The requested diagnostic gas-giant class was not used")
    model.var.verbose = 0
    model.gdat.verbose = False
    if args.reinit_at is not None:
        model.var.nsteps_before_reinit = args.reinit_at
    if args.clear_persistent:
        if args.version != "new" or args.nonpersistent:
            parser.error("--clear-persistent requires --version new without --nonpersistent")
        model.clear_press_temp_edd_profile()
    print(json.dumps({"version": args.version, "persistent": args.version == "new" and not (args.nonpersistent or args.clear_persistent),
                      "model_class": type(model).__name__,
                      "mechanism": temporary_mechanism or "photochem_rxns.yaml",
                      "initialization_seconds": time.monotonic() - start,
                      "equilibrium_time": model.var.equilibrium_time,
                      "conv_longdy": model.var.conv_longdy,
                      "conv_longdydt": model.var.conv_longdydt,
                      "conv_hist_factor": model.var.conv_hist_factor}), flush=True)
    snapshot(model, "initialized")
    model.initialize_robust_stepper(model.wrk.usol)
    snapshot(model, "stepper_initialized")
    previous_temperature = model.var.temperature
    previous_pressure = model.wrk.pressure_hydro
    for i in range(1, args.steps + 1):
        model.var.verbose = args.verbose_level if args.verbose_start <= i <= args.verbose_end else 0
        give_up, converged = model.robust_step()
        detailed = args.detail_start <= i <= args.detail_end
        if detailed:
            temperature = model.var.temperature
            pressure = model.wrk.pressure_hydro
            change = temperature - previous_temperature
            index = int(np.argmax(np.abs(change)))
            print(json.dumps({"detail_step": i,
                              "max_temperature_change": float(change[index]),
                              "temperature_change_layer": index,
                              "max_log10_pressure_change": float(np.max(np.abs(
                                  np.log10(pressure / previous_pressure))))}), flush=True)
        previous_temperature = model.var.temperature
        previous_pressure = model.wrk.pressure_hydro
        if i <= 5 or i % args.interval == 0 or give_up or converged:
            snapshot(model, "stepped", i, converged, give_up)
        if i == args.probe_at:
            same_state = model.wrk.usol.copy()
            before_temperature = model.var.temperature.copy()
            before_pressure = model.wrk.pressure_hydro.copy()
            model.prep_atmosphere(same_state)
            once_temperature = model.var.temperature.copy()
            once_pressure = model.wrk.pressure_hydro.copy()
            model.prep_atmosphere(same_state)
            twice_temperature = model.var.temperature.copy()
            twice_pressure = model.wrk.pressure_hydro.copy()
            print(json.dumps({"probe_step": i,
                              "temperature_first_change": float(np.max(np.abs(
                                  once_temperature - before_temperature))),
                              "temperature_second_change": float(np.max(np.abs(
                                  twice_temperature - once_temperature))),
                              "pressure_first_change": float(np.max(np.abs(
                                  once_pressure - before_pressure))),
                              "pressure_second_change": float(np.max(np.abs(
                                  twice_pressure - once_pressure)))}), flush=True)
            with open("photochem_rxns.yaml", encoding="utf-8") as stream:
                mechanism = yaml.safe_load(stream)
            boundaries = []
            for species in mechanism["species"]:
                thermo = species.get("thermo", {})
                for boundary in thermo.get("temperature-ranges", [])[1:-1]:
                    distance = np.min(np.abs(twice_temperature - boundary))
                    boundaries.append((float(distance), species["name"], float(boundary)))
            nearest = sorted(boundaries)[:5]
            layer_302km = int(np.argmin(np.abs(model.var.z / 1e5 - 302)))
            print(json.dumps({"probe_step": i,
                              "layer_near_302km": layer_302km,
                              "temperature_near_302km": float(twice_temperature[layer_302km]),
                              "nearest_thermo_boundaries": nearest}), flush=True)
        if give_up or converged:
            break
    print(json.dumps({"elapsed_seconds": time.monotonic() - start}), flush=True)
    if temporary_mechanism is not None:
        os.unlink(temporary_mechanism)


if __name__ == "__main__":
    main()
