#!/bin/bash

# Parse command line options for the specified version, force flag, and rebuild flag
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -v|--version) specified_version="$2"; shift;;
        -f|--force) force_flag="true";;
        -r|--rebuild) rebuild_flag="true";;
        *) echo "Error: Unknown parameter passed: $1"; exit 1;;
    esac
    shift
done

# Get the current working directory
current_dir=$(pwd)

# Get the current working directory's dirname
dir_name=$(basename "$current_dir")

# Check if the dirname is "pydock3"
if [ "$dir_name" == "pydock3" ]; then
    # Check if the pyproject.toml file exists
    if [ ! -f "pyproject.toml" ]; then
        echo "Error: pyproject.toml file not found. Exiting."
        exit 2
    fi

    # Get the version from pyproject.toml file
    toml_version=$(grep -oP '(?<=version = ")([^"]*)' pyproject.toml)

    # Check if poetry is installed
    if ! command -v poetry &> /dev/null; then
        echo "Error: poetry is not installed or not in the PATH. Exiting."
        exit 4
    fi

    # Check if pip is installed
    if ! command -v pip &> /dev/null; then
        echo "Error: pip is not installed or not in the PATH. Exiting."
        exit 5
    fi

    # Check if pydock3 is installed
    installed_version=$(pip show pydock3 2>/dev/null | grep -oP "Version: \K.*" || echo "")

    # Exit if no version is specified and pydock3 is already installed unless --force
    if [ -z "$specified_version" ] && [ -n "$installed_version" ] && [ -z "$force_flag" ]; then
        echo "A version of pydock3 is already installed (v${installed_version}). Exiting."
        exit 0
    fi

    if [ -z "$specified_version" ]; then
        specified_version=$toml_version
    else
        if [ -n "$rebuild_flag" ] && [ "$specified_version" != "$toml_version" ]; then
            echo "Error: Cannot build the specified version '${specified_version}' since does not match the version in pyproject.toml (${toml_version}). Exiting."
            exit 8
        fi
    fi

    # Get the dist directory path
    dist_dir="dist/"

    # Find the .whl file in the dist directory for the specified version
    if [ -d "$dist_dir" ]; then
        target_whl_file=$dist_dir/$(ls "$dist_dir" | grep "pydock3-${specified_version}.*\.whl")
    else
        target_whl_file=""
    fi

    #
    if [ "$installed_version" == "$specified_version" ] && [ -z "$force_flag" ]; then
        echo "The specified version '${specified_version}' is already installed. Exiting."
        exit 0
    elif [ -n "$installed_version" ]; then
        pip uninstall -y pydock3
    fi

    # Build the package if the target_whl_file is not found and the specified_version matches the toml_version, or if the rebuild_flag is set
    if [ ! -f "$target_whl_file" ]; then
        if [ "$specified_version" != "$toml_version" ]; then 
            echo "No .whl file found for specified version '${specified_version}' and cannot build new one since specified version does not match the version in pyproject.toml (${toml_version}). Exiting."
            exit 7
        fi
        poetry build -f wheel
        target_whl_file=$dist_dir/$(ls "$dist_dir" | grep "pydock3-${specified_version}.*.whl")
    else 
        if [ -n "$rebuild_flag" ]; then
            poetry build -f wheel
            target_whl_file=$dist_dir/$(ls "$dist_dir" | grep "pydock3-${specified_version}.*.whl")
        fi
    fi

    pip install "$target_whl_file"
    exit 0
else
    echo "Error: The current working directory dirname is NOT 'pydock3'. Exiting."
    exit 6
fi

