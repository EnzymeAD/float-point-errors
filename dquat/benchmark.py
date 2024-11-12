#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import pickle
import numpy as np
import re
import multiprocessing
from matplotlib import pyplot as plt
from matplotlib import rcParams
from tqdm import tqdm

TIMEOUT_SECONDS = 1000
NUM_RUNS = 10


def run_command(
    command, description, capture_output=False, output_file=None, verbose=True, env=None, timeout=TIMEOUT_SECONDS
):
    try:
        if capture_output and output_file:
            with open(output_file, "w") as f:
                subprocess.run(
                    command,
                    stdout=f,
                    stderr=subprocess.STDOUT,
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
        cmd = [executable]
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
    cmd = [executable, "--output-path", output_file]
    run_command(
        cmd, f"Running {executable} and collecting output", capture_output=False, output_file=None, verbose=False
    )


def compute_error(reference_output, test_output):
    return compute_error_geometric_average(reference_output, test_output)


def compute_error_geometric_average(reference_output, test_output):
    ref_matrices_dict = extract_all_matrices(reference_output)
    test_matrices_dict = extract_all_matrices(test_output)

    if ref_matrices_dict is None or test_matrices_dict is None:
        print("Could not extract matrices for error computation.")
        return None

    if set(ref_matrices_dict.keys()) != set(test_matrices_dict.keys()):
        print("Mismatch in iteration numbers between reference and test outputs.")
        return None

    all_relative_diffs = []

    ref_norm = {}
    for iteration in ref_matrices_dict:
        ref_matrix = ref_matrices_dict[iteration]
        print(ref_matrix)
        ref_norm[iteration] = np.linalg.norm(ref_matrix, "fro")

    for iteration in ref_matrices_dict:
        ref_matrix = ref_matrices_dict[iteration]
        test_matrix = test_matrices_dict.get(iteration, None)

        if test_matrix is None:
            print(f"Missing matrix for iteration {iteration}. Skipping.")
            continue

        if ref_matrix.shape != test_matrix.shape:
            print(f"Mismatch in matrix shape at iteration {iteration}. Skipping.")
            continue

        norm_ref = ref_norm.get(iteration, None)
        if norm_ref is None or norm_ref == 0:
            print(f"Reference matrix norm is zero or undefined for iteration {iteration}. Skipping.")
            continue

        diff = ref_matrix - test_matrix
        fro_diff = np.linalg.norm(diff, "fro")
        relative_diff = fro_diff / norm_ref
        all_relative_diffs.append(relative_diff)

    if not all_relative_diffs:
        print("No valid relative differences to compute error.")
        return np.nan

    try:
        log_diffs = [np.log1p(rd) for rd in all_relative_diffs if rd >= 0]
        if not log_diffs:
            return np.nan
        geometric_avg = np.expm1(np.mean(log_diffs))
        relative_error = geometric_avg * 100
        return relative_error
    except (ValueError, OverflowError) as e:
        print(f"Error computing geometric average: {e}")
        return np.nan


def extract_all_matrices(output_file):
    matrices_dict = {}
    current_iteration = None

    with open(output_file, "r") as f:
        lines = f.readlines()

    idx = 0
    while idx < len(lines):
        line = lines[idx]
        iter_match = re.match(r"Test\s+(\d+):", line)
        if iter_match:
            current_iteration = int(iter_match.group(1))
            matrices_dict[current_iteration] = []
            idx += 1
            continue

        if current_iteration is not None:
            if "dquat_dexpmap ^ T value:" in line:
                if idx + 4 >= len(lines):
                    print(f"Insufficient lines for matrix at iteration {current_iteration} in {output_file}")
                    break
                matrix = []
                for j in range(1, 4):
                    vec_line = lines[idx + j].strip()
                    vec = []
                    for val_str in vec_line.split():
                        if val_str.lower() in ["nan", "-nan"]:
                            vec.append(np.nan)
                        else:
                            try:
                                vec.append(float(val_str))
                            except ValueError:
                                vec.append(np.nan)
                    if len(vec) != 4:
                        print(f"Invalid matrix row length at iteration {current_iteration} in {output_file}")
                        vec = [np.nan] * 4
                    matrix.append(vec)
                matrix = np.array(matrix)
                matrices_dict[current_iteration] = matrix
                continue
        idx += 1
    if not matrices_dict:
        print(f"No matrices found in {output_file}")
        return None
    return matrices_dict


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
    (line1,) = ax1.plot(budgets, runtimes, marker="o", linestyle="-", label="Optimized Runtimes", color=color_runtime)
    ax1.tick_params(axis="y", labelcolor=color_runtime)
    ax1.axhline(y=original_runtime, color="blue", linestyle="--", label="Original Runtime")

    ax2 = ax1.twinx()
    color_error = "tab:green"
    ax2.set_ylabel("Relative Errors (%)", color=color_error)
    (line3,) = ax2.plot(
        budgets, adjusted_errors, marker="s", linestyle="-", label="Optimized Relative Errors", color=color_error
    )
    ax2.tick_params(axis="y", labelcolor=color_error)
    ax2.axhline(y=original_error, color="green", linestyle="--", label="Original Relative Error")
    ax2.set_yscale("symlog", linthresh=1e-14)
    ax2.set_ylim(bottom=0)

    ax1.set_title("Computation Cost Budget vs Runtime and Relative Error for dquat Benchmark")
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

    os.makedirs(plots_dir, exist_ok=True)

    plot_filename1 = os.path.join(plots_dir, "runtime_error_plot.png")
    plt.savefig(plot_filename1, bbox_inches="tight", dpi=300)
    plt.close(fig1)
    print(f"First plot saved to {plot_filename1}")

    fig2, ax3 = plt.subplots(figsize=(10, 8))

    ax3.set_xlabel("Runtimes (seconds)")
    ax3.set_ylabel("Relative Errors (%)")
    ax3.set_title(f"Pareto Front of Optimized dquat Programs ({prefix})")

    scatter1 = ax3.scatter(runtimes, adjusted_errors, label="Optimized Programs", color="blue")

    if original_runtime is not None and original_error is not None:
        scatter2 = ax3.scatter(
            original_runtime,
            original_error,
            marker="x",
            color="red",
            s=100,
            label="Original Program",
        )

    points = np.array(list(zip(runtimes, adjusted_errors)))
    sorted_indices = np.argsort(points[:, 0])
    sorted_points = points[sorted_indices]

    pareto_front = [sorted_points[0]]
    for point in sorted_points[1:]:
        if point[1] <= pareto_front[-1][1]:
            pareto_front.append(point)

    pareto_front = np.array(pareto_front)

    (line_pareto,) = ax3.plot(
        pareto_front[:, 0], pareto_front[:, 1], linestyle="-", color="purple", label="Pareto Front"
    )
    ax3.set_yscale("log")

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
        ncol=len(pareto_lines),
        borderaxespad=0.0,
        frameon=False,
    )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)

    plot_filename2 = os.path.join(plots_dir, f"pareto_front_plot_{prefix}.{output_format}")
    plt.savefig(plot_filename2, bbox_inches="tight", dpi=300)
    plt.close(fig2)
    print(f"Second plot saved to {plot_filename2}")


def compute_geometric_average(differences):
    valid_diffs = [1 + d for d in differences if d >= 0]
    if not valid_diffs:
        return np.nan
    log_diffs = [np.log1p(d) for d in valid_diffs]
    return np.expm1(np.mean(log_diffs))


def accuracy_task(executable, reference_output="dquat_gold.txt"):
    m = re.match(r".*dquat-fpopt-(.*)\.exe", executable)
    if m:
        budget = int(m.group(1))
    else:
        print(f"Could not extract budget from executable name: {executable}")
        return None, None
    output_file = os.path.join("outputs", f"output_budget_{budget}.txt")
    collect_output(executable, output_file)
    error = compute_error(reference_output, output_file)
    if error is not None:
        print(f"Relative error for budget {budget}: {error}%")
    else:
        print(f"Relative error for budget {budget} could not be computed.")

    os.remove(output_file)
    return budget, error


def runtime_task(executable, reference_output="dquat_gold.txt"):
    m = re.match(r".*dquat-fpopt-(.*)\.exe", executable)
    if m:
        budget = int(m.group(1))
    else:
        print(f"Could not extract budget from executable name: {executable}")
        return None, None
    runtime = measure_runtime(executable, num_runs=NUM_RUNS)
    if runtime is not None:
        print(f"Measured runtime for budget {budget}: {runtime}")
    else:
        print(f"Runtime for budget {budget} could not be measured.")
    return budget, runtime


def main():
    parser = argparse.ArgumentParser(description="Measure dquat experiments and plot results")
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
    if not os.path.exists(tmp_dir):
        print(f"Temporary directory '{tmp_dir}' does not exist.")
        sys.exit(1)
    executables = [f for f in os.listdir(tmp_dir) if f.startswith("dquat-fpopt-") and f.endswith(".exe")]
    budgets = []
    executable_paths = []
    for exe in executables:
        m = re.match(r"dquat-fpopt-(.*)\.exe", exe)
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

    original_executable = "./dquat.exe"
    original_output_file = os.path.join("outputs", "output_original.txt")
    reference_output = "dquat_gold.txt"

    if not os.path.exists(reference_output):
        print(f"Reference output file '{reference_output}' does not exist.")
        sys.exit(1)

    if not os.path.exists(original_output_file):
        collect_output(original_executable, original_output_file)

    original_error = compute_error(reference_output, original_output_file)
    if original_error is not None:
        print(f"Relative error for the original binary: {original_error}%")
    else:
        print("Relative error for the original binary could not be computed.")

    original_runtime = measure_runtime(original_executable, num_runs=NUM_RUNS)

    print("Starting accuracy measurements...")
    errors = {}
    with multiprocessing.Pool(processes=32) as pool:
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
    plot_results(data["budgets"], data["runtimes"], data["errors"], data["original_error"], data["original_runtime"])


if __name__ == "__main__":
    main()
