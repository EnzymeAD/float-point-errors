#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import pickle
import numpy as np
import numpy.linalg as la
import re
from statistics import mean
import multiprocessing
from matplotlib import pyplot as plt
from matplotlib import rcParams
from tqdm import tqdm
import math

NUM_RUNS = 1
SIZE = 27000
ITER = 932
TIMEOUT_SECONDS = 1000


def run_command(
    command, description, capture_output=False, output_file=None, verbose=True, env=None, timeout=TIMEOUT_SECONDS
):
    try:
        if capture_output and output_file:
            with open(output_file, "w") as f:
                result = subprocess.run(
                    command,
                    stdout=f,
                    text=True,
                    check=True,
                    env=env,
                    timeout=timeout,
                )
                return None
        elif capture_output:
            result = subprocess.run(command, capture_output=True, text=True, check=True, env=env, timeout=timeout)
            return result.stdout
        else:
            if verbose:
                subprocess.check_call(command, env=env, timeout=timeout)
            else:
                subprocess.check_call(
                    command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env, timeout=timeout
                )
    except subprocess.TimeoutExpired:
        print(f"Timeout after {timeout} seconds during: {description}")
        return None
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
        cmd = [executable, "-i", str(ITER)]
        result = run_command(cmd, f"Running {executable} (run {i+1}/{num_runs})", capture_output=True, verbose=False)
        if result is None:
            runtimes.append(np.nan)
            continue
        runtime_line = next((line for line in result.split("\n") if "Elapsed time" in line), None)
        if runtime_line:
            match = re.search(r"Elapsed time\s*=\s*(\S+)", runtime_line)
            if match:
                try:
                    runtime = float(match.group(1))
                    runtimes.append(runtime)
                except ValueError:
                    print(f"Invalid runtime value extracted from output of {executable}")
                    runtimes.append(np.nan)
            else:
                print(f"Could not parse runtime from output of {executable}")
                runtimes.append(np.nan)
        else:
            print(f"Could not find 'Elapsed time' line in output of {executable}")
            runtimes.append(np.nan)
    if not any(~np.isnan(runtimes)):
        print(f"No valid runtimes collected for {executable}.")
        return np.nan
    average_runtime = np.nanmean(runtimes)
    print(f"Average runtime for {executable}: {average_runtime:.6f} seconds")
    return average_runtime


def collect_output(executable, output_file):
    cmd = [executable, "-p", "-i", str(ITER), "-e"]
    run_command(
        cmd, f"Running {executable} and collecting output", capture_output=True, output_file=output_file, verbose=False
    )


def compute_error(reference_output, test_output, iteration):
    ref_energy = extract_energy(reference_output, iteration)
    print(f"=== Measuring error for {test_output} ===")
    test_energy = extract_energy(test_output, iteration)

    if ref_energy is None or test_energy is None:
        print("Could not extract energy values for error computation.")
        return None

    if ref_energy == 0:
        print("Reference energy value is zero, cannot compute relative error.")
        return None

    relative_error = (abs(ref_energy - test_energy) / abs(ref_energy)) * 100
    return relative_error


def extract_energy(output_file, iteration):
    with open(output_file, "r") as f:
        lines = f.readlines()

    iter_count = 0
    for idx, line in enumerate(lines):
        if line.startswith("cycle = "):
            iter_count += 1
            if iter_count == iteration:
                if idx + 1 < len(lines):
                    energy_line = lines[idx + 1].strip()
                    if "Energy (size" in energy_line:
                        m = re.search(r"Energy \(size (\d+)\): \[(.*)\]", energy_line)
                        if m:
                            energy_values_str = m.group(2)
                            first_val_str = energy_values_str.split(",")[0].strip()

                            if first_val_str.lower() in ["nan", "-nan"]:
                                energy_value = np.nan
                            else:
                                try:
                                    energy_value = float(first_val_str)
                                except ValueError:
                                    energy_value = np.nan

                            return energy_value
                        else:
                            print("Regex did not match for energy line.")
                    else:
                        print("Energy line does not contain 'Energy (size'.")
                else:
                    print(f"Reached end of file without finding energy matrix for iteration {iteration}.")
                return np.nan
    print(f"Iteration {iteration} not found in the output file.")
    return np.nan


def plot_results(
    budgets,
    runtimes,
    errors,
    original_error,
    original_runtime,
    prefix="optimized",
    output_format="png",
    plots_dir="plots",
):
    rcParams["font.family"] = "serif"
    rcParams["font.serif"] = ["Linux Libertine"]
    rcParams["font.size"] = 20
    rcParams["axes.titlesize"] = 24
    rcParams["axes.labelsize"] = 20
    rcParams["xtick.labelsize"] = 18
    rcParams["ytick.labelsize"] = 18
    rcParams["legend.fontsize"] = 18

    budgets = np.array(budgets)
    runtimes = np.array(runtimes)
    errors = np.array(errors)

    valid_indices = ~np.isnan(runtimes)
    budgets = budgets[valid_indices]
    runtimes = runtimes[valid_indices]
    errors = errors[valid_indices]

    max_error = 100
    adjusted_errors = np.copy(errors)
    adjusted_errors[np.isnan(adjusted_errors)] = max_error
    adjusted_errors[adjusted_errors > max_error] = max_error

    assert np.all(adjusted_errors >= 0), "Relative errors contain negative values."

    fig1, ax1 = plt.subplots(figsize=(10, 8))

    color_runtime = "tab:blue"
    ax1.set_xlabel("Computation Cost Budget")
    ax1.set_ylabel("Runtimes (seconds)", color=color_runtime)
    scatter1 = ax1.scatter(budgets, runtimes, marker="o", label="Optimized Runtimes", color=color_runtime, s=5)
    ax1.tick_params(axis="y", labelcolor=color_runtime)
    ax1.axhline(y=original_runtime, color="blue", linestyle="--", label="Original Runtime")

    ax2 = ax1.twinx()
    color_error = "tab:green"
    ax2.set_ylabel("Relative Errors (%)", color=color_error)
    scatter2 = ax2.scatter(
        budgets, adjusted_errors, marker="s", label="Optimized Relative Errors", color=color_error, s=5
    )
    ax2.tick_params(axis="y", labelcolor=color_error)
    ax2.axhline(y=original_error, color="green", linestyle="--", label="Original Relative Error")
    ax2.axhline(y=0, color="gold", linestyle="-", label="Zero Relative Error")
    ax2.set_yscale("symlog", linthresh=1e-14)
    ax2.set_ylim(bottom=-1e-14)

    ax1.grid(True)

    lines = [scatter1, scatter2]
    labels = [line.get_label() for line in lines]
    lines.append(plt.Line2D([0], [0], color="blue", linestyle="--"))
    labels.append("Original Runtime")
    lines.append(plt.Line2D([0], [0], color="green", linestyle="--"))
    labels.append("Original Relative Error")
    lines.append(plt.Line2D([0], [0], color="gold", linestyle="-"))

    ax1.legend(
        lines,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=2,
        borderaxespad=0.0,
        frameon=False,
    )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)

    os.makedirs(plots_dir, exist_ok=True)

    plot_filename1 = os.path.join(plots_dir, f"lulesh_runtime_error.{output_format}")
    plt.savefig(plot_filename1, bbox_inches="tight", dpi=300)
    plt.close(fig1)
    print(f"First plot saved to {plot_filename1}")

    fig2, ax3 = plt.subplots(figsize=(10, 8))

    ax3.set_xlabel("Runtimes (seconds)")
    ax3.set_ylabel("Relative Errors (%)")
    scatter1 = ax3.scatter(runtimes, adjusted_errors, label="Optimized Programs", color="blue", s=1)

    if original_runtime is not None and original_error is not None:
        scatter2 = ax3.scatter(
            original_runtime,
            original_error,
            marker="x",
            color="red",
            s=50,
            label="Original Program",
        )

    points = np.array(list(zip(runtimes, adjusted_errors)))
    sorted_indices = np.argsort(points[:, 0])
    sorted_points = points[sorted_indices]

    pareto_front = [sorted_points[0]]
    for point in sorted_points[1:]:
        if point[1] < pareto_front[-1][1]:
            pareto_front.append(point)

    pareto_front = np.array(pareto_front)

    (line_pareto,) = ax3.plot(
        pareto_front[:, 0], pareto_front[:, 1], linestyle="-", color="purple", label="Pareto Front", linewidth=1
    )
    ax3.set_yscale("symlog", linthresh=1e-14)
    ax3.set_ylim(bottom=-1e-14)

    ax3.grid(True)

    pareto_lines = [scatter1, line_pareto]
    pareto_labels = [scatter1.get_label(), line_pareto.get_label()]
    if original_runtime is not None and original_error is not None:
        pareto_lines.append(scatter2)
        pareto_labels.append(scatter2.get_label())


    ax3.legend(
        pareto_lines,
        pareto_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=2,
        borderaxespad=0.0,
        frameon=False,
    )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)

    plot_filename2 = os.path.join(plots_dir, f"lulesh_pareto.{output_format}")
    plt.savefig(plot_filename2, bbox_inches="tight", dpi=300)
    plt.close(fig2)
    print(f"Second plot saved to {plot_filename2}")


def compute_geometric_average(differences):
    valid_diffs = [1 + d for d in differences if d >= 0]
    if not valid_diffs:
        return np.nan
    log_diffs = [np.log1p(d) for d in valid_diffs]
    return np.expm1(np.mean(log_diffs))


def analyze_data(data, thresholds=None):
    budgets = data["budgets"]
    runtimes = data["runtimes"]
    errors = data["errors"]
    original_runtime = data["original_runtime"]
    original_error = data["original_error"]
    print("Original relative error: ", original_error)

    if thresholds is None:
        thresholds = [
            0,
            1e-16,
            1e-15,
            1e-14,
            1e-12,
            1e-10,
            1e-9,
            1e-8,
            1e-6,
            1e-4,
            1e-2,
            1e-1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.9,
            1,
        ]

    if original_error == 0:
        original_digits = 16
    elif original_error is not None:
        original_digits = min(-math.log10(original_error / 100), 16)
    else:
        original_digits = None

    print(f"\nOriginal program has {original_digits:.2f} decimal digits of accuracy.")

    digits_list = []
    for err in errors:
        if err == 0:
            digits = 16
        elif err is not None:
            digits = min(-math.log10(err / 100), 16)
        else:
            digits = None
        digits_list.append(digits)

    accuracy_improvements = []
    for digits in digits_list:
        if digits is not None and original_digits is not None:
            improvement = digits - original_digits
            accuracy_improvements.append(improvement)
        else:
            accuracy_improvements.append(None)

    max_improvement = None
    for improvement in accuracy_improvements:
        if improvement is not None and improvement > 0:
            if max_improvement is None or improvement > max_improvement:
                max_improvement = improvement
    if max_improvement is None:
        max_improvement = 0.0

    print(f"Maximum accuracy improvement: {max_improvement:.2f} decimal digits")

    min_runtime_ratios = {}
    for threshold in thresholds:
        min_ratio = None
        for err, runtime in zip(errors, runtimes):
            if err is not None and runtime is not None and err <= threshold * 100:
                if original_runtime == 0:
                    continue
                runtime_ratio = runtime / original_runtime
                if min_ratio is None or runtime_ratio < min_ratio:
                    min_ratio = runtime_ratio
        if min_ratio is not None:
            min_runtime_ratios[threshold] = min_ratio

    overall_runtime_improvements = {}
    for threshold in thresholds:
        ratio = min_runtime_ratios.get(threshold)
        if ratio is not None:
            percentage_improvement = (1 - ratio) * 100
            overall_runtime_improvements[threshold] = percentage_improvement
        else:
            overall_runtime_improvements[threshold] = None

    print("\nPercentage of runtime improvements while allowing some level of relative error:")
    for threshold in thresholds:
        percentage_improvement = overall_runtime_improvements[threshold]
        if percentage_improvement is not None:
            print(
                f"Allowed relative error ≤ {threshold}: {percentage_improvement:.2f}% runtime reduction / {1 / (1 - percentage_improvement / 100):.2f}x speedup"
            )
        else:
            print(f"Allowed relative error ≤ {threshold}: No data")


def accuracy_task(executable, reference_output="lulesh_gold.txt"):
    m = re.match(r".*ser-single-forward-fpopt-(.*)\.exe", executable)
    if m:
        budget = int(m.group(1))
    else:
        print(f"Could not extract budget from executable name: {executable}")
        return None, None
    output_file = os.path.join("outputs", f"output_budget_{budget}.txt")
    collect_output(executable, output_file)
    error = compute_error(reference_output, output_file, iteration=ITER)
    if error is not None:
        print(f"Relative error for budget {budget}: {error}%")
    else:
        print(f"Relative error for budget {budget} could not be computed.")

    os.remove(output_file)
    return budget, error


def runtime_task(executable, reference_output="lulesh_gold.txt"):
    m = re.match(r".*ser-single-forward-fpopt-(.*)\.exe", executable)
    if m:
        budget = int(m.group(1))
    else:
        print(f"Could not extract budget from executable name: {executable}")
        return None, None
    runtime = measure_runtime(executable, NUM_RUNS)
    if runtime is not None:
        print(f"Measured runtime for budget {budget}: {runtime}")
    else:
        print(f"Runtime for budget {budget} could not be measured.")
    return budget, runtime


def main():
    parser = argparse.ArgumentParser(description="Measure LULESH experiments and plot results")
    parser.add_argument("--plot-only", action="store_true", help="Only plot using data from saved file")
    parser.add_argument(
        "--output-format", type=str, choices=["png", "pdf"], default="png", help="Output format for plots"
    )
    args = parser.parse_args()
    plot_only = args.plot_only
    output_format = args.output_format

    if plot_only:
        if not os.path.exists("measurements.pkl"):
            print("Measurements file 'measurements.pkl' does not exist.")
            sys.exit(1)
        with open("measurements.pkl", "rb") as f:
            data = pickle.load(f)
        plot_results(
            data["budgets"],
            data["runtimes"],
            data["errors"],
            data["original_error"],
            data["original_runtime"],
            output_format=output_format,
        )
        analyze_data(data)
        sys.exit(0)

    if not os.path.exists("plots"):
        os.makedirs("plots")
    if not os.path.exists("outputs"):
        os.makedirs("outputs")

    tmp_dir = "tmp"
    if not os.path.exists(tmp_dir):
        print(f"Temporary directory '{tmp_dir}' does not exist.")
        sys.exit(1)
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
    if not budgets:
        print("No valid executables found.")
        sys.exit(1)
    budgets, executable_paths = zip(*sorted(zip(budgets, executable_paths)))

    original_executable = "./ser-single-forward.exe"
    original_output_file = os.path.join("outputs", "output_original.txt")
    reference_output = "lulesh_gold.txt"

    if not os.path.exists(reference_output):
        print(f"Reference output file '{reference_output}' does not exist.")
        sys.exit(1)

    if not os.path.exists(original_output_file):
        collect_output(original_executable, original_output_file)

    original_error = compute_error(reference_output, original_output_file, iteration=ITER)
    if original_error is not None:
        print(f"Relative error for the original binary: {original_error}%")
    else:
        print("Relative error for the original binary could not be computed.")

    original_runtime = measure_runtime(original_executable, NUM_RUNS)

    print("Starting accuracy measurements...")
    errors = {}
    cpu_count = min(160, multiprocessing.cpu_count())
    with multiprocessing.Pool(processes=cpu_count) as pool:
        results = pool.map(accuracy_task, executable_paths)
        for result in results:
            if result and result[1] is not None:
                budget = result[0]
                errors[budget] = result[1]

    print("Starting runtime measurements...")
    runtimes = {}
    for exe, budget in tqdm(zip(executable_paths, budgets), total=len(budgets), desc="Measuring runtimes"):
        result = runtime_task(exe)
        if result and result[1] is not None:
            runtimes[budget] = result[1]

    data = {
        "budgets": budgets,
        "errors": [errors.get(budget, np.nan) for budget in budgets],
        "runtimes": [runtimes.get(budget, np.nan) for budget in budgets],
        "original_error": original_error if original_error is not None else np.nan,
        "original_runtime": original_runtime if original_runtime is not None else np.nan,
    }

    with open("measurements.pkl", "wb") as f:
        pickle.dump(data, f)
    print("Measurements saved to 'measurements.pkl'")
    plot_results(
        data["budgets"],
        data["runtimes"],
        data["errors"],
        data["original_error"],
        data["original_runtime"],
        output_format=output_format,
    )
    analyze_data(data)


if __name__ == "__main__":
    main()
