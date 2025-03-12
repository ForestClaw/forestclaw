#! /usr/bin/env python3

import os
import subprocess
import shutil
import time
import argparse

def run_ctest(build_dir, rerun_failed=False):
    """Run ctest in the specified build directory."""
    if rerun_failed:
        result = subprocess.run(['ctest', '--rerun-failed'], cwd=build_dir)
    else:
        result = subprocess.run(['ctest'], cwd=build_dir)
    return result.returncode == 0

def handle_actual_files(src_dir):
    """Recursively rename .actual files by removing the .actual extension."""
    actual_files_found = False
    for root, _, files in os.walk(src_dir):
        for file in files:
            if file.endswith('.actual'):
                actual_files_found = True
                actual_file_path = os.path.join(root, file)
                original_file_path = os.path.join(root, file[:-7])  # Remove the .actual extension
                shutil.move(actual_file_path, original_file_path)
    return actual_files_found

def main(build_dir, src_dir):
    """Run ctest repeatedly until no .actual files are produced."""
    first_run = True
    while True:
        print("Running ctest...")
        if run_ctest(build_dir, rerun_failed=not first_run):
            print("All tests passed.")
            break
        print("Handling .actual files...")
        if not handle_actual_files(src_dir):
            print("No .actual files found.")
            break
        print(".actual files handled. Re-running failed tests...")
        first_run = False
        time.sleep(1)  # Small delay to avoid rapid re-running

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ctest and handle .actual files.')
    parser.add_argument('--build', type=str, required=True, help='Path to the build directory')
    parser.add_argument('--src', type=str, required=True, help='Path to the source directory')
    args = parser.parse_args()

    main(args.build, args.src)