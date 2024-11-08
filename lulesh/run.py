#!/usr/bin/env python3

import os
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from tqdm import tqdm

HOME = "/home/sbrantq"

np.random.seed(42)

def run_command(command, description, log_file, env=None):
    try:
        with open(log_file, "w") as f:
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=True,
                env=env
            )
            f.write(result.stdout)
            return result
    except subprocess.CalledProcessError as e:
        print(f"Error during: {description}")
        print(f"Check the log file: {log_file} for details.")
        with open(log_file, "a") as f:
            f.write(f"\nError: {e}\n")
        return None


def compile_shared_objects(env):
    shared_sources = ["lulesh-viz.cc", "lulesh-util.cc", "lulesh-init.cc"]
    shared_objs = []
    for src in shared_sources:
        obj = src.replace('.cc', '.o')
        shared_objs.append(obj)
        cmd_compile_shared = [CXX, "-DOMP_MERGE=0", "-DUSE_MPI=0"] + CXXFLAGS + ["-c", src, "-o", obj]
        log_file = os.path.join("logs", f"compile_shared_{src.replace('.cc', '')}.log")
        run_command(
            cmd_compile_shared,
            f"Compiling shared source file {src}",
            log_file=log_file,
            env=env,
        )
    print("Shared object files compiled successfully.")
    return shared_objs


def compile_binary(budget, env, tmp_dir):
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
        log_file=log_file_compile,
        env=env,
    )

    if compile_output is None:
        print(f"Skipping budget {budget} due to compilation failure.")
        return None

    return obj_file


def link_binary(budget, obj_file, shared_objs, env, tmp_dir):
    """Link the compiled object file with shared objects to create the executable."""
    executable = os.path.join(tmp_dir, f"ser-single-forward-fpopt-{budget}.exe")
    cmd_link = (
        [CXX]
        + CXXFLAGS
        + [obj_file] + shared_objs
        + ["-lm", "-o", executable]
    )

    log_file_link = os.path.join("logs", f"link_budget_{budget}.log")
    result = run_command(
        cmd_link,
        f"Linking binary for budget {budget}",
        log_file=log_file_link,
        env=env,
    )

    if result is None:
        print(f"Linking failed for budget {budget}.")
        return None

    return executable


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
        "-I" + os.path.join(HOME, "include"),
        "-L" + os.path.join(HOME, "lib"),
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
        # "-mllvm",
        # "-fpopt-early-prune",
        "-mllvm",
        "-herbie-timeout=500",
        "-mllvm",
        "-herbie-num-threads=16",
        # "-mllvm",
        # "-herbie-disable-taylor",
        "-mllvm",
        "-fpopt-cost-model-path=cm.csv",
    ]

    run_command(["make", "clean"], "Cleaning previous builds", log_file=os.path.join("logs", "make_clean.log"))

    shared_objs = compile_shared_objects(env)

    NUM_TESTED_COSTS = 128
    budget_range = [-570000000000, -380000000000]
    budget_lower = min(budget_range)
    budget_upper = max(budget_range)

    budgets = np.random.uniform(budget_lower, budget_upper, NUM_TESTED_COSTS).astype(int).tolist()
    budgets.sort()
    print("Testing the following budgets:", budgets)

    # Phase 1: Compile all binaries in parallel
    max_workers = 64
    compiled_objs = []
    failed_budgets = []

    print("=== Phase 1: Compiling Object Files ===")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_budget = {executor.submit(compile_binary, budget, env, "tmp"): budget for budget in budgets}

        with tqdm(total=NUM_TESTED_COSTS, desc="Compiling Binaries") as pbar:
            for future in as_completed(future_to_budget):
                budget = future_to_budget[future]
                try:
                    obj_file = future.result()
                    if obj_file:
                        compiled_objs.append((budget, obj_file))
                    else:
                        failed_budgets.append(budget)
                except Exception as e:
                    print(f"An error occurred during compilation for budget {budget}: {e}")
                    failed_budgets.append(budget)
                finally:
                    pbar.update(1)

    print(f"Compilation phase completed. {len(compiled_objs)} succeeded, {len(failed_budgets)} failed.")

    # Phase 2: Link binaries sequentially
    print("=== Phase 2: Linking Binaries ===")
    with tqdm(total=len(compiled_objs), desc="Linking Binaries") as pbar:
        for budget, obj_file in compiled_objs:
            executable = link_binary(budget, obj_file, shared_objs, env, "tmp")
            if executable:
                print(f"Linked executable for budget {budget}: {executable}")
            else:
                print(f"Linking failed for budget {budget}.")
            pbar.update(1)

    print("All binaries compiled and linked. Executables are saved in the 'tmp' folder.")
    if failed_budgets:
        print("Some budgets failed during compilation or linking:")
        for b in failed_budgets:
            print(f"  - Budget {b}")


if __name__ == "__main__":
    main()
