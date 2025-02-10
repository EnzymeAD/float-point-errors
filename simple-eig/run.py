#!/usr/bin/env python3

import os
import subprocess
import sys
import shutil
import time
import re
import argparse
import matplotlib
import matplotlib.pyplot as plt
import math
import random
import numpy as np
from statistics import mean
import pickle
import seaborn as sns

from tqdm import tqdm, trange
from matplotlib import rcParams
from matplotlib.legend_handler import HandlerTuple
from concurrent.futures import ProcessPoolExecutor, as_completed

PREFIX = "simple-eig"

HOME = "/home/sbrantq"
ENZYME_PATH = os.path.join(HOME, "sync/Enzyme/build-debug/Enzyme/ClangEnzyme-16.so")
LLVM_PATH = os.path.join(HOME, "llvms/llvm16/build/bin")
CXX = os.path.join(LLVM_PATH, "clang++")

CXXFLAGS = [
    "-O3",
    "-I" + os.path.join(HOME, "include"),
    "-L" + os.path.join(HOME, "lib"),
    "-I/usr/include/c++/11",
    "-I/usr/include/x86_64-linux-gnu/c++/11",
    "-L/usr/lib/gcc/x86_64-linux-gnu/11",
    "-fno-exceptions",
    f"-fpass-plugin={ENZYME_PATH}",
    "-Xclang",
    "-load",
    "-Xclang",
    ENZYME_PATH,
    "-lmpfr",
    "-ffast-math",
    "-fno-finite-math-only",
    "-fuse-ld=lld",
]

FPOPTFLAGS_BASE = [
    "-mllvm",
    "--enzyme-enable-fpopt",
    "-mllvm",
    "--enzyme-print-herbie",
    "-mllvm",
    "--enzyme-print-fpopt",
    "-mllvm",
    "--fpopt-log-path=example.txt",
    "-mllvm",
    "--fpopt-target-func-regex=example",
    "-mllvm",
    "--fpopt-enable-solver",
    "-mllvm",
    "--fpopt-enable-pt",
    "-mllvm",
    "--fpopt-comp-cost-budget=0",
    "-mllvm",
    "--herbie-num-threads=8",
    "-mllvm",
    "--herbie-timeout=1000",
    "-mllvm",
    "--fpopt-num-samples=1024",
    "-mllvm",
    "--fpopt-cost-model-path=/home/sbrantq/sync/FPBench/microbm/cm.csv",
    "-mllvm",
    "-fpopt-cache-path=cache",
]

SRC = "example.c"
LOGGER = "fp-logger.cpp"
EXE = ["example.exe", "example-logged.exe", "example-fpopt.exe"]
NUM_RUNS = 10
DRIVER_NUM_SAMPLES = 10000000
LOG_NUM_SAMPLES = 10000
MAX_TESTED_COSTS = 999


def run_command(command, description, capture_output=False, output_file=None, verbose=True, timeout=None):
    print(f"=== {description} ===")
    print("Running:", " ".join(command))
    try:
        if capture_output and output_file:
            with open(output_file, "w") as f:
                subprocess.check_call(command, stdout=f, stderr=subprocess.STDOUT, timeout=timeout)
        elif capture_output:
            result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=timeout)
            return result.stdout
        else:
            if verbose:
                subprocess.check_call(command, timeout=timeout)
            else:
                subprocess.check_call(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout)
    except subprocess.TimeoutExpired:
        print(f"Command '{' '.join(command)}' timed out after {timeout} seconds.")
        return
    except subprocess.CalledProcessError as e:
        print(f"Error during: {description}")
        if capture_output and output_file:
            print(f"Check the output file: {output_file} for details.")
        else:
            print(e)
        sys.exit(e.returncode)


def clean(tmp_dir, logs_dir, plots_dir):
    print("=== Cleaning up generated files ===")
    directories = [tmp_dir, logs_dir, plots_dir]
    for directory in directories:
        if os.path.exists(directory):
            shutil.rmtree(directory)
            print(f"Removed directory: {directory}")


def generate_example_cpp(tmp_dir, prefix):
    script = "fpopt-original-driver-generator.py"
    print(f"=== Running {script} ===")
    src_prefixed = os.path.join(tmp_dir, f"{prefix}{SRC}")
    dest_prefixed = os.path.join(tmp_dir, f"{prefix}example.cpp")
    run_command(
        ["python3", script, src_prefixed, dest_prefixed, "example", str(DRIVER_NUM_SAMPLES)],
        f"Generating example.cpp from {SRC}",
    )
    if not os.path.exists(dest_prefixed):
        print(f"Failed to generate {dest_prefixed}.")
        sys.exit(1)
    print(f"Generated {dest_prefixed} successfully.")


def generate_example_logged_cpp(tmp_dir, prefix):
    script = "fpopt-logged-driver-generator.py"
    print(f"=== Running {script} ===")
    src_prefixed = os.path.join(tmp_dir, f"{prefix}{SRC}")
    dest_prefixed = os.path.join(tmp_dir, f"{prefix}example-logged.cpp")
    run_command(
        ["python3", script, src_prefixed, dest_prefixed, "example", str(LOG_NUM_SAMPLES)],
        f"Generating example-logged.cpp from {SRC}",
    )
    if not os.path.exists(dest_prefixed):
        print(f"Failed to generate {dest_prefixed}.")
        sys.exit(1)
    print(f"Generated {dest_prefixed} successfully.")


def generate_example_baseline_cpp(tmp_dir, prefix):
    script = "fpopt-baseline-generator.py"
    print(f"=== Running {script} ===")
    src_prefixed = os.path.join(tmp_dir, f"{prefix}{SRC}")
    dest_prefixed = os.path.join(tmp_dir, f"{prefix}example-baseline.cpp")
    run_command(
        ["python3", script, src_prefixed, dest_prefixed, "example", str(DRIVER_NUM_SAMPLES)],
        f"Generating example-baseline.cpp from {SRC}",
    )
    if not os.path.exists(dest_prefixed):
        print(f"Failed to generate {dest_prefixed}.")
        sys.exit(1)
    print(f"Generated {dest_prefixed} successfully.")


def compile_example_exe(tmp_dir, prefix):
    source = os.path.join(tmp_dir, f"{prefix}example.cpp")
    output = os.path.join(tmp_dir, f"{prefix}example.exe")
    cmd = [CXX, source] + CXXFLAGS + ["-o", output]
    run_command(cmd, f"Compiling {output}")


def compile_example_logged_exe(tmp_dir, prefix):
    source = os.path.join(tmp_dir, f"{prefix}example-logged.cpp")
    output = os.path.join(tmp_dir, f"{prefix}example-logged.exe")
    cmd = [CXX, source, LOGGER] + CXXFLAGS + ["-o", output]
    run_command(cmd, f"Compiling {output}")


def compile_example_baseline_exe(tmp_dir, prefix):
    source = os.path.join(tmp_dir, f"{prefix}example-baseline.cpp")
    output = os.path.join(tmp_dir, f"{prefix}example-baseline.exe")
    cmd = [CXX, source] + CXXFLAGS + ["-o", output]
    run_command(cmd, f"Compiling {output}")


def generate_example_txt(tmp_dir, prefix):
    exe = os.path.join(tmp_dir, f"{prefix}example-logged.exe")
    output = os.path.join(tmp_dir, f"{prefix}example.txt")
    if not os.path.exists(exe):
        print(f"Executable {exe} not found. Cannot generate {output}.")
        sys.exit(1)
    with open(output, "w") as f:
        print(f"=== Running {exe} to generate {output} ===")
        try:
            subprocess.check_call([exe], stdout=f)
        except subprocess.TimeoutExpired:
            print(f"Execution of {exe} timed out.")
            if os.path.exists(exe):
                os.remove(exe)
                print(f"Removed executable {exe} due to timeout.")
            return
        except subprocess.CalledProcessError as e:
            print(f"Error running {exe}")
            sys.exit(e.returncode)


def compile_example_fpopt_exe(tmp_dir, prefix, fpoptflags, output="example-fpopt.exe", verbose=True):
    source = os.path.join(tmp_dir, f"{prefix}example.cpp")
    output_path = os.path.join(tmp_dir, f"{prefix}{output}")
    cmd = [CXX, source] + CXXFLAGS + fpoptflags + ["-o", output_path]
    log_path = os.path.join("logs", f"{prefix}compile_fpopt.log")
    if output == "example-fpopt.exe":
        run_command(
            cmd,
            f"Compiling {output_path} with FPOPTFLAGS",
            capture_output=True,
            output_file=log_path,
            verbose=verbose,
        )
    else:
        run_command(
            cmd,
            f"Compiling {output_path} with FPOPTFLAGS",
            verbose=verbose,
        )


def parse_critical_comp_costs(tmp_dir, prefix, log_path="compile_fpopt.log"):
    print(f"=== Parsing critical computation costs from {log_path} ===")
    full_log_path = os.path.join("logs", f"{prefix}{log_path}")
    if not os.path.exists(full_log_path):
        print(f"Log file {full_log_path} does not exist.")
        sys.exit(1)
    with open(full_log_path, "r") as f:
        content = f.read()

    pattern = r"\*\*\* Critical Computation Costs \*\*\*(.*?)\*\*\* End of Critical Computation Costs \*\*\*"
    match = re.search(pattern, content, re.DOTALL)
    if not match:
        print("Critical Computation Costs block not found in the log.")
        sys.exit(1)

    costs_str = match.group(1).strip()
    costs = [int(cost) for cost in costs_str.split(",") if re.fullmatch(r"-?\d+", cost.strip())]
    print(f"Parsed computation costs: {costs}")

    if not costs:
        print("No valid computation costs found to sample.")
        sys.exit(1)

    num_to_sample = min(MAX_TESTED_COSTS, len(costs))

    sampled_costs = random.sample(costs, num_to_sample)

    sampled_costs_sorted = sorted(sampled_costs)

    print(f"Sampled computation costs (sorted): {sampled_costs_sorted}")

    return sampled_costs_sorted


def measure_runtime(tmp_dir, prefix, executable, num_runs=NUM_RUNS):
    print(f"=== Measuring runtime for {executable} ===")
    runtimes = []
    exe_path = os.path.join(tmp_dir, f"{prefix}{executable}")
    for i in trange(1, num_runs + 1):
        try:
            result = subprocess.run([exe_path], capture_output=True, text=True, check=True, timeout=300)
            output = result.stdout
            match = re.search(r"Total runtime: ([\d\.]+) seconds", output)
            if match:
                runtime = float(match.group(1))
                runtimes.append(runtime)
            else:
                print(f"Could not parse runtime from output on run {i}")
                sys.exit(1)
        except subprocess.TimeoutExpired:
            print(f"Execution of {exe_path} timed out on run {i}")
            if os.path.exists(exe_path):
                os.remove(exe_path)
                print(f"Removed executable {exe_path} due to timeout.")
            return None
        except subprocess.CalledProcessError as e:
            print(f"Error running {exe_path} on run {i}")
            sys.exit(e.returncode)
    if runtimes:
        average_runtime = mean(runtimes)
        print(f"Average runtime for {prefix}{executable}: {average_runtime:.6f} seconds")
        return average_runtime
    else:
        print(f"No successful runs for {prefix}{executable}")
        return None


def get_values_file_path(tmp_dir, prefix, binary_name):
    return os.path.join(tmp_dir, f"{prefix}{binary_name}-values.txt")


def generate_example_values(tmp_dir, prefix):
    binary_name = "example.exe"
    exe = os.path.join(tmp_dir, f"{prefix}{binary_name}")
    output_values_file = get_values_file_path(tmp_dir, prefix, binary_name)
    cmd = [exe, "--output-path", output_values_file]
    run_command(cmd, f"Generating function values from {binary_name}", verbose=False, timeout=300)


def generate_values(tmp_dir, prefix, binary_name):
    exe = os.path.join(tmp_dir, f"{prefix}{binary_name}")
    values_file = get_values_file_path(tmp_dir, prefix, binary_name)
    cmd = [exe, "--output-path", values_file]
    run_command(cmd, f"Generating function values from {binary_name}", verbose=False, timeout=300)


def compile_golden_exe(tmp_dir, prefix):
    source = os.path.join(tmp_dir, f"{prefix}golden.cpp")
    output = os.path.join(tmp_dir, f"{prefix}golden.exe")
    cmd = [CXX, source] + CXXFLAGS + ["-o", output]
    run_command(cmd, f"Compiling {output}")


def generate_golden_values(tmp_dir, prefix):
    script = "fpopt-golden-driver-generator.py"
    src_prefixed = os.path.join(tmp_dir, f"{prefix}{SRC}")
    dest_prefixed = os.path.join(tmp_dir, f"{prefix}golden.cpp")
    cur_prec = 128
    max_prec = 4096
    PREC_step = 128
    prev_output = None
    output_values_file = get_values_file_path(tmp_dir, prefix, "golden.exe")
    while cur_prec <= max_prec:
        run_command(
            ["python3", script, src_prefixed, dest_prefixed, str(cur_prec), "example", str(DRIVER_NUM_SAMPLES)],
            f"Generating golden.cpp with PREC={cur_prec}",
        )
        if not os.path.exists(dest_prefixed):
            print(f"Failed to generate {dest_prefixed}.")
            sys.exit(1)
        print(f"Generated {dest_prefixed} successfully.")

        compile_golden_exe(tmp_dir, prefix)

        exe = os.path.join(tmp_dir, f"{prefix}golden.exe")
        cmd = [exe, "--output-path", output_values_file]
        run_command(cmd, f"Generating golden values with PREC={cur_prec}", verbose=False)

        if not os.path.exists(output_values_file):
            print(f"Failed to generate golden values at PREC={cur_prec} due to timeout.")
            return  # Assume golden values do not exist

        with open(output_values_file, "r") as f:
            output = f.read()

        if output == prev_output:
            print(f"Golden values converged at PREC={cur_prec}")
            break
        else:
            prev_output = output
            cur_prec += PREC_step
    else:
        print(f"Failed to converge golden values up to PREC={max_prec}")
        sys.exit(1)


def get_avg_rel_error(tmp_dir, prefix, golden_values_file, binaries):
    with open(golden_values_file, "r") as f:
        golden_values = [float(line.strip()) for line in f]

    errors = {}
    for binary in binaries:
        values_file = get_values_file_path(tmp_dir, prefix, binary)
        if not os.path.exists(values_file):
            print(f"Values file {values_file} does not exist. Skipping error calculation for {binary}.")
            errors[binary] = None
            continue
        with open(values_file, "r") as f:
            values = [float(line.strip()) for line in f]
        if len(values) != len(golden_values):
            print(f"Number of values in {values_file} does not match golden values")
            sys.exit(1)

        valid_errors = []
        for v, g in zip(values, golden_values):
            if math.isnan(v) or math.isnan(g):
                continue
            if g == 0:
                continue
            error = abs((v - g) / g) * 100
            valid_errors.append(error)

        if not valid_errors:
            print(f"No valid data to compute rel error for binary {binary}. Setting rel error to None.")
            errors[binary] = None
            continue

        try:
            log_sum = sum(math.log1p(e) for e in valid_errors)
            geo_mean = math.expm1(log_sum / len(valid_errors))
            errors[binary] = geo_mean
        except OverflowError:
            print(
                f"Overflow error encountered while computing geometric mean for binary {binary}. Setting rel error to None."
            )
            errors[binary] = None
        except ZeroDivisionError:
            print(f"No valid errors to compute geometric mean for binary {binary}. Setting rel error to None.")
            errors[binary] = None

    return errors


def plot_results(
    plots_dir, prefix, budgets, runtimes, errors, original_runtime=None, original_error=None, output_format="pdf"
):
    print(f"=== Plotting results to {output_format.upper()} file ===")

    # Filter out invalid data
    data = list(zip(budgets, runtimes, errors))
    # Removing budget 0 since we know it is the original program for this benchmark
    filtered_data = [(b, r, e) for b, r, e in data if b != 0 and r is not None and e is not None]

    if not filtered_data:
        print("No valid data to plot.")
        return

    budgets, runtimes, errors = zip(*filtered_data)

    # Print all runtime–relative error pairs sorted by budget
    print(f"Original Program: {original_runtime:.6f} seconds, {original_error:.6e}%")
    print("Runtime and Relative Error pairs (sorted by Computation Cost Budget):")
    for b, r, e in zip(budgets, runtimes, errors):
        print(f"  Budget: {b}, Runtime: {r:.6f} seconds, Relative Error: {e:.6e}%")

    # -------------------------------------------------------------------------
    # Use Seaborn's "whitegrid" style/palette
    # -------------------------------------------------------------------------
    print(plt.style.available)
    sns.set_theme(style="whitegrid")

    # (Optional) Adjust some font sizes
    rcParams["font.family"] = "serif"
    rcParams["font.serif"] = ["Linux Libertine"]
    # rcParams["font.size"] = 24
    # rcParams["axes.titlesize"] = 24
    rcParams["axes.labelsize"] = 24
    rcParams["xtick.labelsize"] = 20
    rcParams["ytick.labelsize"] = 20
    rcParams["legend.fontsize"] = 20

    # -------------------------------------------------------------------------
    # First Plot: Budget vs. (Runtime and Relative Error)
    # -------------------------------------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(10, 8))

    color_runtime = "C0"  # Let Seaborn pick color C0
    color_error = "C1"  # Let Seaborn pick color C1

    # Plot runtimes on left axis
    ax1.set_xlabel("Computation Cost Budget")
    ax1.set_ylabel("Runtimes (seconds)", color=color_runtime)
    (line1,) = ax1.step(
        budgets, runtimes, marker="o", linestyle="-", label="Optimized Runtimes", color=color_runtime, where="post"
    )
    if original_runtime is not None:
        line2 = ax1.axhline(y=original_runtime, color=color_runtime, linestyle="--", label="Original Runtime")
    ax1.tick_params(axis="y", labelcolor=color_runtime)

    # Plot errors on right axis
    ax2 = ax1.twinx()
    ax2.set_ylabel("Relative Errors (%)", color=color_error)
    (line3,) = ax2.step(
        budgets,
        errors,
        marker="s",
        linestyle="-.",
        label="Optimized Relative Errors",
        color=color_error,
        where="post",
    )
    if original_error is not None:
        line4 = ax2.axhline(y=original_error, color=color_error, linestyle="--", label="Original Relative Error")
    ax2.tick_params(axis="y", labelcolor=color_error)
    # Use a "symlog" scale with a tiny linthresh to handle near-zero errors
    ax2.set_yscale("symlog", linthresh=1e-14)
    ax2.set_ylim(bottom=0)

    ax1.set_title(f"Computation Cost Budget vs Runtime\nand Relative Error ({prefix[:-1]})")
    ax1.grid(True, linestyle=":", alpha=0.7)

    # Build combined legend for the first plot
    lines = [line1, line3]
    labels = [line1.get_label(), line3.get_label()]
    if original_runtime is not None:
        lines.append(line2)
        labels.append(line2.get_label())
    if original_error is not None:
        lines.append(line4)
        labels.append(line4.get_label())

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

    plot_filename1 = os.path.join(plots_dir, f"runtime_plot_{prefix[:-1]}.{output_format}")
    plt.savefig(plot_filename1, bbox_inches="tight", dpi=300)
    plt.close(fig1)
    print(f"First plot saved to {plot_filename1}")

    # -------------------------------------------------------------------------
    # Second Plot: Pareto Front (Runtimes vs. Relative Errors)
    # -------------------------------------------------------------------------
    fig2, ax3 = plt.subplots(figsize=(10, 8))

    blue = "#2287E6"
    yellow = "#FFBD59"
    red = "#FF6666"
    big_star_size = 300
    small_star_size = 250
    triangle_size = 200

    ax3.set_xlabel("Runtimes (seconds)")
    ax3.set_ylabel("Relative Errors (%)")
    # ax3.set_title(f"Pareto Front of Optimized Programs ({prefix[:-1]})")
    ax3.set_yscale("log")
    ax3.grid(True, linestyle=":", alpha=0.7)

    # -------------------------------------------------------------------------
    # Identify Pareto-optimal points (including the original program if provided)
    # -------------------------------------------------------------------------
    # First, build the array of points from optimized programs.
    optimization_points = np.array(list(zip(runtimes, errors)))
    # If an original program is provided, add it to the set of points.
    if original_runtime is not None and original_error is not None:
        all_points = np.vstack([optimization_points, [original_runtime, original_error]])
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
    line_pareto = None
    if len(pareto_points) > 0:
        sorted_idx = np.argsort(pareto_points[:, 0])
        pareto_line_points = pareto_points[sorted_idx]
        (line_pareto,) = ax3.step(
            pareto_line_points[:, 0],
            pareto_line_points[:, 1],
            linestyle="--",
            color=blue,
            label="Pareto Front",
            where="post",
            zorder=-1,
        )

    # -------------------------------------------------------------------------
    # Markers for special vs. other budgets
    # -------------------------------------------------------------------------
    special_budgets = {-30590000, 10450000, -29580000}
    special_handle = None
    other_handle = None
    for b, rt, err in zip(budgets, runtimes, errors):
        if b in special_budgets:
            h = ax3.scatter(rt, err, marker="*", s=big_star_size, color=blue, zorder=10)
            if special_handle is None:
                special_handle = h
        else:
            h = ax3.scatter(rt, err, marker="*", s=small_star_size, color=yellow, zorder=5)
            if other_handle is None:
                other_handle = h

    # Plot the original program as a separate marker (if provided)
    if original_runtime is not None and original_error is not None:
        scatter_original = ax3.scatter(
            original_runtime,
            original_error,
            marker="^",
            color=red,
            s=triangle_size,
            label="Original Program",
            zorder=10,
        )

    # -------------------------------------------------------------------------
    # Build legend: Combine the two optimized markers into one entry
    # -------------------------------------------------------------------------
    legend_elements = []
    legend_labels = []

    # Combine special and other handles into a tuple for a single legend entry.
    if special_handle is not None or other_handle is not None:
        if special_handle is not None and other_handle is not None:
            optimized_handle = (special_handle, other_handle)
        else:
            optimized_handle = special_handle if special_handle is not None else other_handle
        legend_elements.append(optimized_handle)
        legend_labels.append("Optimized")

    if original_runtime is not None and original_error is not None:
        legend_elements.append(scatter_original)
        legend_labels.append("Original")
    if line_pareto is not None:
        legend_elements.append(line_pareto)
        legend_labels.append("Pareto Front")

    ax3.legend(
        legend_elements,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=len(legend_elements),
        borderaxespad=0.0,
        frameon=False,
        handler_map={tuple: HandlerTuple(ndivide=None)},
    )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)

    plot_filename2 = os.path.join(plots_dir, f"pareto_front_plot_{prefix[:-1]}.{output_format}")
    plt.savefig(plot_filename2, bbox_inches="tight", dpi=300)
    plt.close(fig2)
    print(f"Second plot saved to {plot_filename2}")


def build_all(tmp_dir, logs_dir, prefix):
    generate_example_cpp(tmp_dir, prefix)
    generate_example_logged_cpp(tmp_dir, prefix)
    # generate_example_baseline_cpp(tmp_dir, prefix)
    compile_example_exe(tmp_dir, prefix)
    compile_example_logged_exe(tmp_dir, prefix)
    # compile_example_baseline_exe(tmp_dir, prefix)
    generate_example_txt(tmp_dir, prefix)
    fpoptflags = []
    for flag in FPOPTFLAGS_BASE:
        if flag.startswith("--fpopt-log-path="):
            fpoptflags.append(f"--fpopt-log-path=tmp/{prefix}example.txt")
        else:
            fpoptflags.append(flag)
    compile_example_fpopt_exe(tmp_dir, prefix, fpoptflags, output="example-fpopt.exe")
    print("=== Initial build process completed successfully ===")


def measure_baseline_runtime(tmp_dir, prefix, num_runs=NUM_RUNS):
    executable = f"example-baseline.exe"
    avg_runtime = measure_runtime(tmp_dir, prefix, executable, num_runs)
    return avg_runtime


def process_cost(args):
    cost, tmp_dir, prefix = args

    print(f"\n=== Processing computation cost budget: {cost} ===")
    fpoptflags = []
    for flag in FPOPTFLAGS_BASE:
        if flag.startswith("--fpopt-comp-cost-budget="):
            fpoptflags.append(f"--fpopt-comp-cost-budget={cost}")
        elif flag.startswith("--fpopt-log-path="):
            fpoptflags.append(f"--fpopt-log-path=tmp/{prefix}example.txt")
        else:
            fpoptflags.append(flag)

    output_binary = f"example-fpopt-{cost}.exe"

    compile_example_fpopt_exe(tmp_dir, prefix, fpoptflags, output=output_binary, verbose=False)

    generate_values(tmp_dir, prefix, output_binary)

    return cost, output_binary


def benchmark(tmp_dir, logs_dir, prefix, plots_dir, num_parallel=1):
    costs = parse_critical_comp_costs(tmp_dir, prefix)

    original_avg_runtime = measure_runtime(tmp_dir, prefix, "example.exe", NUM_RUNS)
    original_runtime = original_avg_runtime

    if original_runtime is None:
        print("Original binary timed out. Proceeding as if it doesn't exist.")
        return

    generate_example_values(tmp_dir, prefix)

    generate_golden_values(tmp_dir, prefix)

    golden_values_file = get_values_file_path(tmp_dir, prefix, "golden.exe")
    example_binary = "example.exe"
    rel_errs_example = get_avg_rel_error(tmp_dir, prefix, golden_values_file, [example_binary])
    rel_err_example = rel_errs_example[example_binary]
    print(f"Average Rel Error for {prefix}example.exe: {rel_err_example}")

    data_tuples = []

    args_list = [(cost, tmp_dir, prefix) for cost in costs]

    if num_parallel == 1:
        for args in args_list:
            cost, output_binary = process_cost(args)
            data_tuples.append((cost, output_binary))
    else:
        with ProcessPoolExecutor(max_workers=num_parallel) as executor:
            future_to_cost = {executor.submit(process_cost, args): args[0] for args in args_list}
            for future in as_completed(future_to_cost):
                cost = future_to_cost[future]
                try:
                    cost_result, output_binary = future.result()
                    data_tuples.append((cost_result, output_binary))
                except Exception as exc:
                    print(f"Cost {cost} generated an exception: {exc}")

    data_tuples_sorted = sorted(data_tuples, key=lambda x: x[0])
    sorted_budgets, sorted_optimized_binaries = zip(*data_tuples_sorted) if data_tuples_sorted else ([], [])

    # Measure runtimes serially based on sorted budgets
    sorted_runtimes = []
    for cost, output_binary in zip(sorted_budgets, sorted_optimized_binaries):
        avg_runtime = measure_runtime(tmp_dir, prefix, output_binary, NUM_RUNS)
        if avg_runtime is not None:
            sorted_runtimes.append(avg_runtime)
        else:
            print(f"Skipping cost {cost} due to runtime measurement failure.")
            sorted_runtimes.append(None)

    errors_dict = get_avg_rel_error(tmp_dir, prefix, golden_values_file, sorted_optimized_binaries)
    sorted_errors = []
    for binary in sorted_optimized_binaries:
        sorted_errors.append(errors_dict.get(binary))
        print(f"Average rel error for {binary}: {errors_dict.get(binary)}")

    sorted_budgets = list(sorted_budgets)
    sorted_runtimes = list(sorted_runtimes)
    sorted_errors = list(sorted_errors)

    data = {
        "budgets": sorted_budgets,
        "runtimes": sorted_runtimes,
        "errors": sorted_errors,
        "original_runtime": original_runtime,
        "original_error": rel_err_example,
    }
    data_file = os.path.join(tmp_dir, f"{prefix}benchmark_data.pkl")
    with open(data_file, "wb") as f:
        pickle.dump(data, f)
    print(f"Benchmark data saved to {data_file}")

    plot_results(
        plots_dir,
        prefix,
        sorted_budgets,
        sorted_runtimes,
        sorted_errors,
        original_runtime=original_runtime,
        original_error=rel_err_example,
    )


def plot_from_data(tmp_dir, plots_dir, prefix, output_format="pdf"):
    data_file = os.path.join(tmp_dir, f"{prefix}benchmark_data.pkl")
    if not os.path.exists(data_file):
        print(f"Data file {data_file} does not exist. Cannot plot.")
        sys.exit(1)
    with open(data_file, "rb") as f:
        data = pickle.load(f)
    plot_results(
        plots_dir,
        prefix,
        data["budgets"],
        data["runtimes"],
        data["errors"],
        original_runtime=data["original_runtime"],
        original_error=data["original_error"],
        output_format=output_format,
    )


def analyze_all_data(tmp_dir, thresholds=None):
    prefixes = []
    data_list = []

    for filename in os.listdir(tmp_dir):
        if filename.endswith("benchmark_data.pkl"):
            data_file = os.path.join(tmp_dir, filename)
            with open(data_file, "rb") as f:
                data = pickle.load(f)
            prefix = filename[: -len("benchmark_data.pkl")]
            prefixes.append(prefix)
            data_list.append((prefix, data))

    print("Number of tested FPBench functions: ", len(data_list))
    if not data_list:
        print("No benchmark data files found in the tmp directory.")
        return

    print(f"Analyzing data for prefixes: {', '.join(prefixes)}\n")

    if thresholds is None:
        thresholds = [0, 1e-10, 1e-9, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.9, 1]

    max_accuracy_improvements = {}  # Per benchmark
    min_runtime_ratios = {threshold: {} for threshold in thresholds}

    original_digits = []
    for prefix, data in data_list:
        budgets = data["budgets"]
        runtimes = data["runtimes"]
        errors = data["errors"]
        original_runtime = data["original_runtime"]
        original_error = data["original_error"]

        if original_error == 0:
            example_digits = 17
            # print(f"Original program of {prefix} has zero relative error! Using 17 digits of accuracy.")
        else:
            example_digits = min(-math.log10(original_error / 100), 17)
            # print(f"Original program of {prefix} has {example_digits:.2f} digits of accuracy.")

        original_digits.append(example_digits)

        digits_list = []
        for err in errors:
            if err == 0:
                digits_list.append(17)
                # print(f"Optimized program of {prefix} has zero relative error! Using 17 digits of accuracy.")
            elif err is not None:
                digits = min(-math.log10(err / 100), 17)
                digits_list.append(digits)
                # print(f"Optimized program of {prefix} has {digits:.2f} digits of accuracy.")
            else:
                digits_list.append(None)

        accuracy_improvements = []
        for digits in digits_list:
            if digits is not None and example_digits is not None:
                improvement = digits - example_digits
                accuracy_improvements.append(improvement)
            else:
                accuracy_improvements.append(None)

        # Find the maximum accuracy improvement for this benchmark
        max_improvement = None
        for improvement in accuracy_improvements:
            if improvement is not None and improvement > 0:
                if max_improvement is None or improvement > max_improvement:
                    max_improvement = improvement

        max_accuracy_improvements[prefix] = max_improvement or 0.0

        # For each threshold, find the minimum runtime ratio for this benchmark
        for threshold in thresholds:
            min_ratio = None
            for err, runtime in zip(errors, runtimes):
                if err is not None and runtime is not None and err <= threshold * 100:
                    # print(f"Threshold: {threshold}, Error: {err}, Runtime: {runtime}")
                    runtime_ratio = runtime / original_runtime
                    if min_ratio is None or runtime_ratio < min_ratio:
                        min_ratio = runtime_ratio
            if min_ratio is not None:
                min_runtime_ratios[threshold][prefix] = min_ratio
                # print(f"Threshold: {threshold}, maximum runtime improvement for {prefix}: {(1 - min_ratio) * 100:.2f}%")

    log_sum = sum(math.log1p(digits) for digits in original_digits)
    geo_mean = math.expm1(log_sum / len(original_digits))
    print(f"Original programs have {geo_mean:.2f} decimal digits of accuracy on average.")
    print(f"Original programs have {max(original_digits):.2f} decimal digits of accuracy at most.")

    # Now compute the geometric mean of minimum runtime ratios per threshold
    overall_runtime_improvements = {}
    for threshold in thresholds:
        ratios = min_runtime_ratios[threshold].values()
        # print(f"\nThreshold: {threshold}, Number of valid runtime ratios: {len(ratios)}")
        if ratios:
            log_sum = sum(math.log1p(min(1, ratio)) for ratio in ratios)
            geo_mean_ratio = math.expm1(log_sum / len(ratios))
            percentage_improvement = (1 - geo_mean_ratio) * 100
            overall_runtime_improvements[threshold] = percentage_improvement
        else:
            overall_runtime_improvements[threshold] = None

    # ratios1 = min_runtime_ratios[1e-06].keys()
    # ratios2 = min_runtime_ratios[1e-04].keys()

    # for i in ratios2:
    #     if i not in ratios1:
    #         print(i)

    # Print maximum accuracy improvements per benchmark
    print("Maximum accuracy improvements (in number of bits) per benchmark:")
    for prefix in prefixes:
        improvement = max_accuracy_improvements.get(prefix)
        if improvement:
            print(f"{prefix}: {improvement:.2f} decimal digits")
        else:
            print(f"{prefix}: No improvement")

    improvements = list(max_accuracy_improvements.values())

    # Compute geometric mean of maximum accuracy improvements
    positive_improvements = [impr for impr in improvements if impr > 0]

    if not positive_improvements:
        print("\nNo positive accuracy improvements available to compute geometric mean.")
    else:
        try:
            log_sum = sum(math.log(impr) for impr in positive_improvements)
            geo_mean = math.exp(log_sum / len(positive_improvements))
            print(f"\nGeometric mean of maximum accuracy improvements: {geo_mean:.2f} decimal digits")
        except ValueError as e:
            print(f"\nError in computing geometric mean: {e}")

    # Print overall runtime improvements per threshold
    print("\nGeometric average percentage of runtime improvements while allowing some level of relative error:")
    for threshold in thresholds:
        percentage_improvement = overall_runtime_improvements[threshold]
        if percentage_improvement is not None:
            print(f"Allowed relative error ≤ {threshold}: {percentage_improvement:.2f}% runtime improvement")
        else:
            print(f"Allowed relative error ≤ {threshold}: No data")


def build_with_benchmark(tmp_dir, logs_dir, plots_dir, prefix, num_parallel=1):
    build_all(tmp_dir, logs_dir, prefix)
    benchmark(tmp_dir, logs_dir, prefix, plots_dir, num_parallel)


def remove_cache_dir():
    cache_dir = "cache"
    if os.path.exists(cache_dir):
        shutil.rmtree(cache_dir)
        print("=== Removed existing cache directory ===")


def copy_source_to_tmp(tmp_dir, prefix):
    """
    Copies 'simple-eig.c' from the current directory into tmp,
    renaming it to '<prefix>example.c' as expected by the build.
    """
    os.makedirs(tmp_dir, exist_ok=True)
    src_file = PREFIX + ".c"
    dst_file = os.path.join(tmp_dir, f"{prefix}{SRC}")
    if not os.path.exists(dst_file):
        try:
            shutil.copy(src_file, dst_file)
            print(f"Copied {src_file} to {dst_file}")
        except FileNotFoundError:
            print(f"Error: {src_file} not found in the current directory.")
            sys.exit(1)
    else:
        print(f"Source file already exists at {dst_file}")


def main():
    parser = argparse.ArgumentParser(description="Run the example C code with prefix handling.")
    parser.add_argument("--clean", action="store_true", help="Clean up generated files")
    parser.add_argument("--build", action="store_true", help="Build all components")
    parser.add_argument("--benchmark", action="store_true", help="Run benchmark")
    parser.add_argument("--all", action="store_true", help="Build and run benchmark")
    parser.add_argument("--plot-only", action="store_true", help="Plot results from existing data")
    parser.add_argument("--output-format", type=str, default="pdf", help="Output format for plots (e.g., png, pdf)")
    parser.add_argument("--analytics", action="store_true", help="Run analytics on saved data")
    parser.add_argument("--disable-preopt", action="store_true", help="Disable Enzyme preoptimization")
    parser.add_argument(
        "--num-parallel", type=int, default=16, help="Number of parallel processes to use (default: 16)"
    )
    args = parser.parse_args()

    global FPOPTFLAGS_BASE
    if args.disable_preopt:
        FPOPTFLAGS_BASE.extend(["-mllvm", "--enzyme-preopt=0"])

    prefix = PREFIX
    if not prefix.endswith("-"):
        prefix += "-"

    tmp_dir = "tmp"
    logs_dir = "logs"
    plots_dir = "plots"

    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    # remove_cache_dir()
    copy_source_to_tmp(tmp_dir, prefix)

    if args.clean:
        clean(tmp_dir, logs_dir, plots_dir)
        sys.exit(0)
    elif args.build:
        build_all(tmp_dir, logs_dir, prefix)
        sys.exit(0)
    elif args.benchmark:
        benchmark(tmp_dir, logs_dir, prefix, plots_dir, num_parallel=args.num_parallel)
        sys.exit(0)
    elif args.plot_only:
        plot_from_data(tmp_dir, plots_dir, prefix, output_format=args.output_format)
        sys.exit(0)
    elif args.analytics:
        analyze_all_data(tmp_dir)
        sys.exit(0)
    elif args.all:
        build_with_benchmark(tmp_dir, logs_dir, plots_dir, prefix, num_parallel=args.num_parallel)
        sys.exit(0)
    else:
        build_with_benchmark(tmp_dir, logs_dir, plots_dir, prefix, num_parallel=args.num_parallel)


if __name__ == "__main__":
    main()
