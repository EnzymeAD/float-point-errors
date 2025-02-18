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
import seaborn as sns
from matplotlib.legend_handler import HandlerTuple
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

    relative_error = abs(ref_energy - test_energy) / abs(ref_energy)
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
    sns.set_theme(style="whitegrid")
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

    valid_indices = ~np.isnan(runtimes) & ~np.isnan(errors) & (runtimes <= 15)
    budgets = budgets[valid_indices]
    runtimes = runtimes[valid_indices]
    errors = errors[valid_indices]

    max_error = 1.0
    adjusted_errors = np.copy(errors)
    adjusted_errors[np.isnan(adjusted_errors)] = max_error
    adjusted_errors[adjusted_errors > max_error] = max_error

    assert np.all(adjusted_errors >= 0), "Relative errors contain negative values."

    fig1, ax1 = plt.subplots(figsize=(10, 8))

    color_runtime = "C0"
    color_error = "C1"

    ax1.set_xlabel("Computation Cost Budget")
    ax1.set_ylabel("Runtimes (seconds)", color=color_runtime)
    (line1,) = ax1.step(
        budgets, runtimes, marker="o", linestyle="-", label="Optimized Runtimes", color=color_runtime, where="post"
    )
    if original_runtime is not None:
        line2 = ax1.axhline(y=original_runtime, color=color_runtime, linestyle="--", label="Original Runtime")
    ax1.tick_params(axis="y", labelcolor=color_runtime)

    ax2 = ax1.twinx()
    ax2.set_ylabel("Relative Errors", color=color_error)
    (line3,) = ax2.step(
        budgets,
        adjusted_errors,
        marker="s",
        linestyle="-.",
        label="Optimized Relative Errors",
        color=color_error,
        where="post",
    )
    if original_error is not None:
        line4 = ax2.axhline(y=original_error, color=color_error, linestyle="--", label="Original Relative Error")
    ax2.tick_params(axis="y", labelcolor=color_error)
    ax2.set_yscale("symlog", linthresh=1e-14)
    ax2.set_ylim(bottom=0)

    ax1.set_title(f"Computation Cost Budget vs Runtime\nand Relative Error ({prefix})")
    ax1.grid(True, linestyle=":", alpha=0.7)

    # Build a combined legend for the first plot
    lines = [line1, line3]
    labels = [line1.get_label(), line3.get_label()]
    if original_runtime is not None:
        lines.append(line2)
        labels.append("Original Runtime")
    if original_error is not None:
        lines.append(line4)
        labels.append("Original Relative Error")
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

    # -------------------------------------------------------------------------
    # Second Plot: Pareto Front (Runtimes vs. Relative Errors) with broken y-axis
    # We remove the region from 1e-11 to 1e0 by using two subplots.
    # -------------------------------------------------------------------------
    # Create two subplots that share the x-axis
    fig2, (ax_upper, ax_lower) = plt.subplots(
        2, 1, sharex=True, figsize=(10, 8), gridspec_kw={"height_ratios": [0.2, 0.8]}
    )

    # Use a symlog scale on both axes
    for ax in (ax_upper, ax_lower):
        ax.set_yscale("symlog", linthresh=1e-14, linscale=0.5)

    # TODO: Fix
    ax_upper.set_ylim(1e-1, 500)
    ax_lower.set_ylim(-1e-15, 1e-11)

    # Remove the spines between the two plots
    ax_upper.spines["bottom"].set_visible(False)
    ax_lower.spines["top"].set_visible(False)
    ax_upper.tick_params(labelbottom=False)  # hide x tick labels on the upper plot

    # Set common labels

    # Define colors and marker sizes
    blue = "#2287E6"
    yellow = "#FFBD59"
    red = "#FF6666"
    big_star_size = 100
    small_star_size = 10
    triangle_size = 100

    # Identify Pareto-optimal points (including the original program if provided)
    optimization_points = np.array(list(zip(runtimes, adjusted_errors, budgets)))
    if original_runtime is not None and original_error is not None:
        all_points = np.vstack([optimization_points, [original_runtime, original_error, 0]])
    else:
        all_points = optimization_points

    n_points = len(all_points)
    pareto_optimal = np.zeros(n_points, dtype=bool)
    current_best_error = float("inf")
    indices_by_runtime = np.argsort(all_points[:, 0])
    for idx in indices_by_runtime:
        if all_points[idx, 1] < current_best_error:
            pareto_optimal[idx] = True
            current_best_error = all_points[idx, 1]

    pareto_points = all_points[pareto_optimal]
    non_pareto_points = all_points[~pareto_optimal]

    print("\n=== Pareto Optimal Points ===")
    for rt, err, b in pareto_points:
        print(f"Budget: {b}, Runtime: {rt:.6f} seconds, Relative Error: {err:.6e}")

    # Draw the Pareto front as a step plot.
    # We plot on each axis so that only the portion in its y-range is visible.
    line_pareto = None
    if len(pareto_points) > 0:
        sorted_idx = np.argsort(pareto_points[:, 0])
        pareto_line_points = pareto_points[sorted_idx]
        # Plot on the upper panel (with label)
        line_pareto = ax_upper.step(
            pareto_line_points[:, 0],
            pareto_line_points[:, 1],
            linestyle="--",
            color=blue,
            label="Pareto Front",
            where="post",
            zorder=-1,
        )[0]
        # Plot on the lower panel without a label
        ax_lower.step(
            pareto_line_points[:, 0],
            pareto_line_points[:, 1],
            linestyle="--",
            color=blue,
            where="post",
            zorder=-1,
        )

    # Mark optimized points.
    # For points that belong to a "special" budget (i.e. those on the Pareto front),
    # we use a larger marker.
    special_budgets = {b for _, _, b in pareto_points}
    special_handle = None
    other_handle = None
    for b, rt, err in zip(budgets, runtimes, adjusted_errors):
        if b in special_budgets:
            # Plot on upper panel (add label only the first time)
            if special_handle is None:
                h = ax_upper.scatter(rt, err, marker="*", s=big_star_size, color=blue, zorder=10, label="Optimized")
                special_handle = h
            else:
                ax_upper.scatter(rt, err, marker="*", s=big_star_size, color=blue, zorder=10)
            # Plot on lower panel (no label)
            ax_lower.scatter(rt, err, marker="*", s=big_star_size, color=blue, zorder=10)
        else:
            if other_handle is None:
                h = ax_upper.scatter(rt, err, marker="*", s=small_star_size, color=yellow, zorder=5, label="Optimized")
                other_handle = h
            else:
                ax_upper.scatter(rt, err, marker="*", s=small_star_size, color=yellow, zorder=5)
            ax_lower.scatter(rt, err, marker="*", s=small_star_size, color=yellow, zorder=5)

    # Plot the original program with a distinct marker (if provided)
    if original_runtime is not None and original_error is not None:
        scatter_original = ax_upper.scatter(
            original_runtime,
            original_error,
            marker="^",
            color=red,
            s=triangle_size,
            label="Original Program",
            zorder=10,
        )
        ax_lower.scatter(
            original_runtime,
            original_error,
            marker="^",
            color=red,
            s=triangle_size,
            zorder=10,
        )

    optimized_handle = special_handle if special_handle is not None else other_handle
    legend_elements = []
    legend_labels = []
    if optimized_handle is not None:
        legend_elements.append(optimized_handle)
        legend_labels.append("Optimized")
    if original_runtime is not None and original_error is not None:
        legend_elements.append(scatter_original)
        legend_labels.append("Original")
    if line_pareto is not None:
        legend_elements.append(line_pareto)
        legend_labels.append("Pareto Front")

    fig2.legend(
        legend_elements,
        legend_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.1),
        ncol=len(legend_elements),
        # borderaxespad=0.0,
        frameon=False,
    )

    ax_lower.set_xlabel("Runtimes (seconds)")
    fig2.text(0.0, 0.6, "Relative Errors", va="center", ha="center", rotation="vertical", fontsize=20)

    # Add break markers (diagonal lines) to indicate the break in the y-axis.
    # d = 0.015  # size of diagonal lines in axes coordinates
    # kwargs = dict(color="k", clip_on=False, lw=1)
    # On the upper axes, draw at the bottom
    # ax_upper.plot((-d, +d), (-d, +d), transform=ax_upper.transAxes, **kwargs)
    # ax_upper.plot((1 - d, 1 + d), (-d, +d), transform=ax_upper.transAxes, **kwargs)
    # # On the lower axes, draw at the top
    # ax_lower.plot((-d, +d), (1 - d, 1 + d), transform=ax_lower.transAxes, **kwargs)
    # ax_lower.plot((1 - d, 1 + d), (1 - d, 1 + d), transform=ax_lower.transAxes, **kwargs)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)

    plot_filename2 = os.path.join(plots_dir, f"lulesh_pareto.{output_format}")
    plt.savefig(plot_filename2, bbox_inches="tight", dpi=300)
    plt.close(fig2)
    print(f"Second plot saved to {plot_filename2}")


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

    min_runtime_ratios = {}
    for threshold in thresholds:
        min_ratio = None
        for err, runtime in zip(errors, runtimes):
            if err is not None and runtime is not None and err <= threshold:
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
        print(f"Relative error for budget {budget}: {error}")
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

    collect_output(original_executable, original_output_file)

    original_error = compute_error(reference_output, original_output_file, iteration=ITER)
    os.remove(original_output_file)
    if original_error is not None:
        print(f"Relative error for the original binary: {original_error}")
    else:
        print("Relative error for the original binary could not be computed.")

    original_runtime = measure_runtime(original_executable, NUM_RUNS)

    print("Starting accuracy measurements...")
    errors = {}
    cpu_count = min(64, multiprocessing.cpu_count())
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
