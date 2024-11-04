#!/usr/bin/env python3

import os
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from tqdm import tqdm

HOME = "/home/sbrantq"

np.random.seed(42)

def run_command(command, description, capture_output=False, output_file=None, verbose=True, env=None):
    if verbose:
        print(f"=== {description} ===")
        print("Running:", " ".join(command))
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
        if verbose:
            print(f"Error during: {description}")
            if capture_output and output_file:
                print(f"Check the output file: {output_file} for details.")
            else:
                print(e)
        return None


def compile_shared_objects(env):
    """Compile shared object files once."""
    shared_sources = ["lulesh-viz.cc", "lulesh-util.cc", "lulesh-init.cc"]
    cmd_compile_shared = [CXX, "-DOMP_MERGE=0", "-DUSE_MPI=0"] + CXXFLAGS + ["-c"] + shared_sources
    run_command(
        cmd_compile_shared,
        "Compiling shared source files",
        capture_output=False,
        verbose=False,
        env=env,
    )
    print("Shared object files compiled successfully.")


def build_optimized_binary(budget, env, tmp_dir):
    """Compile lulesh.cc with a specific FPOpt budget and link with shared objects."""
    description = f"Building optimized LULESH binary with budget {budget}"
    fpopt_flags = FPOPTFLAGS_BASE.copy()
    fpopt_flags += ["-mllvm", f"-fpopt-comp-cost-budget={budget}"]

    obj_file = os.path.join(tmp_dir, f"lulesh-fpopt-{budget}.o")
    cmd_compile_fpopt = (
        [CXX, "-DOMP_MERGE=0", "-DUSE_MPI=0", "-DFORWARD=1"]
        + CXXFLAGS
        + ["-c", "lulesh.cc"]
        + fpopt_flags
        + ["-o", obj_file]
    )

    log_file_compile = os.path.join("logs", f"compile_budget_{budget}.log")
    compile_output = run_command(
        cmd_compile_fpopt,
        f"Compiling lulesh.cc with budget {budget}",
        capture_output=True,
        output_file=log_file_compile,
        verbose=False,
        env=env,
    )

    if compile_output is None:
        print(f"Skipping budget {budget} due to compilation failure.")
        return

    executable = os.path.join(tmp_dir, f"ser-single-forward-fpopt-{budget}.exe")
    cmd_link = (
        [CXX]
        + CXXFLAGS
        + [
            obj_file,
            "lulesh-viz.o",
            "lulesh-util.o",
            "lulesh-init.o",
            "-lm",
            "-o",
            executable,
        ]
    )

    try:
        run_command(
            cmd_link,
            f"Linking binary for budget {budget}",
            capture_output=False,
            verbose=False,
            env=env,
        )
        return executable
    except subprocess.CalledProcessError as e:
        print(f"Linking failed for budget {budget}.")
        print(e)
        return


def main():
    parser = argparse.ArgumentParser(description="Compile LULESH optimized binaries with FPOpt")
    args = parser.parse_args()

    os.makedirs("logs", exist_ok=True)
    os.makedirs("tmp", exist_ok=True)

    env = os.environ.copy()
    env["ENZYME_PATH"] = os.path.join(HOME, "sync/Enzyme/build/Enzyme/ClangEnzyme-15.so")
    env["CLANG_PATH"] = os.path.join(HOME, "llvms/llvm15/build/bin")

    env["OPENMP_PATH"] = os.path.join(env["CLANG_PATH"], "../projects/openmp/runtime/src")
    env["MPI_PATH"] = "/usr/lib/x86_64-linux-gnu/openmpi/lib"
    env["OPENMP_LIB"] = os.path.join(env["CLANG_PATH"], "../lib/libomp.so")

    global CXX, CXXFLAGS, FPOPTFLAGS_BASE

    CXX = os.path.join(env["CLANG_PATH"], "clang++")
    CXXFLAGS = [
        "-O3",
        "-I.",
        "-Wall",
        "-fno-exceptions",
        "-I" + os.path.join(env["HOME"], "include"),
        "-L" + os.path.join(env["HOME"], "lib"),
        "-I/usr/include/c++/13",
        "-I/usr/include/x86_64-linux-gnu/c++/13",
        "-L/usr/lib/gcc/x86_64-linux-gnu/13",
        "-I" + env["OPENMP_PATH"],
        "-L" + env["MPI_PATH"],
        "-mllvm",
        "-enzyme-loose-types",
        "-mllvm",
        "-enzyme-inline",
        "-fplugin=" + env["ENZYME_PATH"],
        "-lmpi",
        "-ffast-math",
        "-fno-finite-math-only",
    ]

    FPOPTFLAGS_BASE = [
        "-mllvm",
        "-enzyme-enable-fpopt",
        "-mllvm",
        "-fpopt-log-path=lulesh.txt",
        "-mllvm",
        "-fpopt-enable-solver",
        "-mllvm",
        "-fpopt-target-func-regex=LagrangeLeapFrog",
        "-mllvm",
        "-fpopt-enable-pt",
        "-mllvm",
        "-fpopt-num-samples=1000",
        "-mllvm",
        "-fpopt-cost-dom-thres=0.0",
        "-mllvm",
        "-fpopt-acc-dom-thres=0.0",
        "-mllvm",
        "-fpopt-early-prune",
        "-mllvm",
        "-herbie-timeout=120",
        "-mllvm",
        "-herbie-disable-taylor",
        "-mllvm",
        "-fpopt-cost-model-path=cm.csv",
    ]

    run_command(["make", "clean"], "Cleaning previous builds", verbose=False)

    compile_shared_objects(env)

    NUM_TESTED_COSTS = 128
    budget_range = [-500000000000, -200000000000]
    budget_lower = min(budget_range)
    budget_upper = max(budget_range)

    budgets = np.random.uniform(budget_lower, budget_upper, NUM_TESTED_COSTS).astype(int).tolist()
    budgets.sort()
    print("Testing the following budgets:", budgets)

    max_workers = 64
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_budget = {executor.submit(build_optimized_binary, budget, env, "tmp"): budget for budget in budgets}

        with tqdm(total=NUM_TESTED_COSTS, desc="Compiling Binaries") as pbar:
            for future in as_completed(future_to_budget):
                budget = future_to_budget[future]
                try:
                    result = future.result()
                    if result:
                        pass
                except Exception as e:
                    print(f"An error occurred for budget {budget}: {e}")
                finally:
                    pbar.update(1)

    print("All binaries compiled and saved in 'tmp' folder.")


if __name__ == "__main__":
    main()
