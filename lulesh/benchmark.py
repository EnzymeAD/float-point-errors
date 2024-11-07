#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import pickle
import numpy as np
import re
from statistics import mean
from concurrent.futures import ProcessPoolExecutor, as_completed
from matplotlib import pyplot as plt
from matplotlib import rcParams

NUM_RUNS = 5


def run_command(command, description, capture_output=False, output_file=None, verbose=True, env=None):
    try:
        if capture_output and output_file:
            with open(output_file, "w") as f:
                result = subprocess.run(
                    command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True, env=env
                )
                f.write(result.stdout)
                return result.stdout
        elif capture_output:
            result = subprocess.run(command, capture_output=True, text=True, check=True, env=env)
            return result.stdout
        else:
            if verbose:
                subprocess.check_call(command, env=env)
            else:
                subprocess.check_call(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    except subprocess.CalledProcessError as e:
        print(f"Error during: {description}")
        if capture_output and output_file:
            print(f"Check the output file: {output_file} for details.")
        else:
            print(e)
        return None


def measure_runtime(executable, num_runs):
    print(f"=== Measuring runtime for {executable} ===")
    runtimes = []
    for i in range(num_runs):
        cmd = [executable, "-i", "30"]
        result = run_command(cmd, f"Running {executable} (run {i+1}/{num_runs})", capture_output=True, verbose=False)
        if result is None:
            continue
        runtime_line = next((line for line in result.split("\n") if "Elapsed time" in line), None)
        if runtime_line:
            match = re.search(r"Elapsed time\s*=\s*(\S+)", runtime_line)
            if match:
                runtime = float(match.group(1))
                runtimes.append(runtime)
            else:
                print(f"Could not parse runtime from output of {executable}")
        else:
            print(f"Could not find 'Elapsed time' line in output of {executable}")
    if not runtimes:
        print(f"No runtimes collected for {executable}.")
        return None
    average_runtime = mean(runtimes)
    print(f"Average runtime for {executable}: {average_runtime:.6f} seconds")
    return average_runtime


def collect_output(executable, output_file):
    cmd = [executable, "-p", "-i", "30", "-e"]
    result = run_command(
        cmd, f"Running {executable} and collecting output", capture_output=True, output_file=output_file, verbose=False
    )
    return output_file


def compute_error(reference_output, test_output, iteration):
    ref_energy_value = extract_energy_value(reference_output, iteration)
    print(f"=== Measuring error for {test_output} ===")
    test_energy_value = extract_energy_value(test_output, iteration)
    print(f"Reference energy value: {ref_energy_value}")
    print(f"Test energy value: {test_energy_value}")
    if ref_energy_value is None or test_energy_value is None:
        print("Could not extract energy values for error computation.")
        return None
    diff = ref_energy_value - test_energy_value
    if ref_energy_value == 0:
        print("Reference energy value is zero, cannot compute relative error.")
        return None
    relative_error = (abs(diff) / abs(ref_energy_value)) * 100
    return relative_error


def extract_energy_value(output_file, iteration):
    with open(output_file, "r") as f:
        lines = f.readlines()
    energy_value = None
    iter_count = 0
    line_idx = 0
    while line_idx < len(lines):
        line = lines[line_idx]
        if line.startswith("cycle = "):
            iter_count += 1
            if iter_count == iteration:
                energy_line_index = line_idx + 1
                if energy_line_index < len(lines):
                    energy_line = lines[energy_line_index]
                    if energy_line.startswith("Energy (size"):
                        m = re.search(r"Energy \(size \d+\): \[(.*)\]", energy_line)
                        if m:
                            energy_values_str = m.group(1)
                            energy_values = energy_values_str.strip().split(",")
                            if energy_values:
                                energy_value = float(energy_values[0])
                                return energy_value
        line_idx += 1
    return None


def plot_results(budgets, runtimes, errors, original_error, original_runtime):
    rcParams["font.size"] = 20
    rcParams["axes.titlesize"] = 24
    rcParams["axes.labelsize"] = 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    rcParams["legend.fontsize"] = 18

    valid_indices = [i for i in range(len(errors)) if errors[i] is not None and runtimes[i] is not None]
    budgets = [budgets[i] for i in valid_indices]
    runtimes = [runtimes[i] for i in valid_indices]
    errors = [errors[i] for i in valid_indices]

    fig1, ax1 = plt.subplots(figsize=(10, 8))

    color_runtime = "tab:blue"
    ax1.set_xlabel("Computation Cost Budget")
    ax1.set_ylabel("Runtimes (seconds)", color=color_runtime)
    (line1,) = ax1.plot(budgets, runtimes, marker="o", linestyle="-", label="Optimized Runtimes", color=color_runtime)
    ax1.tick_params(axis="y", labelcolor=color_runtime)
    ax1.axhline(y=original_runtime, color="blue", linestyle="--", label="Original Runtime")

    ax2 = ax1.twinx()
    color_error = "tab:green"
    ax2.set_ylabel("Relative Errors (%)", color=color_error)
    (line3,) = ax2.plot(
        budgets, errors, marker="s", linestyle="-", label="Optimized Relative Errors", color=color_error
    )
    ax2.tick_params(axis="y", labelcolor=color_error)
    ax2.axhline(y=original_error, color="green", linestyle="--", label="Original Relative Error")
    ax2.set_yscale("log")
    ax1.set_title("Computation Cost Budget vs Runtime and Relative Error")
    ax1.grid(True)
    lines = [line1, line3]
    labels = [line.get_label() for line in lines]
    lines.append(plt.Line2D([0], [0], color="blue", linestyle="--"))
    labels.append("Original Runtime")
    lines.append(plt.Line2D([0], [0], color="green", linestyle="--"))
    labels.append("Original Relative Error")
    ax1.legend(
        lines,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=3,
        borderaxespad=0.0,
        frameon=False,
    )
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    plot_filename1 = os.path.join("plots", f"runtime_error_plot_lulesh.png")
    plt.savefig(plot_filename1, bbox_inches="tight", dpi=300)
    plt.close(fig1)
    print(f"Plot saved to {plot_filename1}")


def identify_pareto_optimal(costs, errors):
    num_points = len(costs)
    is_pareto = np.ones(num_points, dtype=bool)
    for i in range(num_points):
        if is_pareto[i]:
            is_pareto[is_pareto] = ~(
                (costs[is_pareto] <= costs[i])
                & (errors[is_pareto] <= errors[i])
                & ((costs[is_pareto] < costs[i]) | (errors[is_pareto] < errors[i]))
            )
            is_pareto[i] = True
    return np.where(is_pareto)[0]


def accuracy_task(executable, reference_output="lulesh_gold.txt"):
    m = re.match(r".*ser-single-forward-fpopt-(.*)\.exe", executable)
    if m:
        budget = int(m.group(1))
    else:
        print(f"Could not extract budget from executable name: {executable}")
        return None, None
    output_file = os.path.join("outputs", f"output_budget_{budget}.txt")
    collect_output(executable, output_file)
    error = compute_error(reference_output, output_file, iteration=30)
    print(f"Relative error for budget {budget}: {error}%")
    return budget, error


def runtime_task(executable, reference_output="lulesh_gold.txt"):
    m = re.match(r".*ser-single-forward-fpopt-(.*)\.exe", executable)
    if m:
        budget = int(m.group(1))
    else:
        print(f"Could not extract budget from executable name: {executable}")
        return None, None
    runtime = measure_runtime(executable, NUM_RUNS)
    print(f"Measured runtime for budget {budget}: {runtime}")
    return budget, runtime


def main():
    parser = argparse.ArgumentParser(description="Measure LULESH experiments and plot results")
    parser.add_argument("--plot-only", action="store_true", help="Only plot using data from saved file")
    args = parser.parse_args()
    plot_only = args.plot_only
    if plot_only:
        if not os.path.exists("measurements.pkl"):
            print("Measurements file 'measurements.pkl' does not exist.")
            sys.exit(1)
        with open("measurements.pkl", "rb") as f:
            data = pickle.load(f)
        plot_results(
            data["budgets"], data["runtimes"], data["errors"], data["original_error"], data["original_runtime"]
        )
        sys.exit(0)
    if not os.path.exists("plots"):
        os.makedirs("plots")
    if not os.path.exists("outputs"):
        os.makedirs("outputs")

    tmp_dir = "tmp"
    executables = [f for f in os.listdir(tmp_dir) if f.startswith("ser-single-forward-fpopt-") and f.endswith(".exe")]
    budgets = []
    executable_paths = []
    for exe in executables:
        m = re.match(r"ser-single-forward-fpopt-(.*)\.exe", exe)
        if m:
            budget_str = m.group(1)
            try:
                budget = int(budget_str)
                budgets.append(budget)
                executable_paths.append(os.path.join(tmp_dir, exe))
            except ValueError:
                continue
    budgets, executable_paths = zip(*sorted(zip(budgets, executable_paths)))

    original_executable = "./ser-single-forward.exe"
    original_output_file = os.path.join("outputs", "output_original.txt")
    reference_output = "lulesh_gold.txt"

    if not os.path.exists(reference_output):
        print(f"Reference output file '{reference_output}' does not exist.")
        sys.exit(1)

    if not os.path.exists(original_output_file):
        run_command(
            [original_executable, "-p", "-i", "30", "-e"],
            "Running original binary and collecting output",
            capture_output=True,
            output_file=original_output_file,
            verbose=False,
        )

    original_error = compute_error(reference_output, original_output_file, iteration=30)
    print(f"Relative error for the original binary: {original_error}%")

    original_runtime = measure_runtime(original_executable, NUM_RUNS)

    print("Starting accuracy measurements...")
    errors = {}
    with ProcessPoolExecutor(64) as executor:
        futures = {executor.submit(accuracy_task, exe): budget for exe, budget in zip(executable_paths, budgets)}
        for future in as_completed(futures):
            budget = futures[future]
            result = future.result()
            if result:
                errors[budget] = result[1]

    print("Starting runtime measurements...")
    runtimes = {}
    for exe, budget in zip(executable_paths, budgets):
        result = runtime_task(exe)
        if result:
            runtimes[budget] = result[1]

    data = {
        "budgets": budgets,
        "errors": [errors.get(budget) for budget in budgets],
        "runtimes": [runtimes.get(budget) for budget in budgets],
        "original_error": original_error,
        "original_runtime": original_runtime,
    }

    with open("measurements.pkl", "wb") as f:
        pickle.dump(data, f)
    print("Measurements saved to 'measurements.pkl'")
    plot_results(data["budgets"], data["runtimes"], data["errors"], data["original_error"], data["original_runtime"])


if __name__ == "__main__":
    main()
