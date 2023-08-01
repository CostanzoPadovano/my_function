#!/bin/bash

# Function to create the analysis folder structure
create_analysis_folder() {
    # Create necessary directories for the analysis project structure
    # The -p flag ensures that parent directories are created if they don't exist

    # Data directory: Holds input data for the analysis
    mkdir -p "$1/data/raw"        # Raw data (if applicable)
    mkdir -p "$1/data/processed"  # Processed and cleaned data

    # Analysis notebooks and scripts
    mkdir -p "$1/notebooks"      # Jupyter notebooks or other analysis scripts
    mkdir -p "$1/scripts"        # Supporting scripts or utility functions

    # Reports and visualizations
    mkdir -p "$1/reports"        # Analysis reports, visualizations, and presentations

    # Results directory: For output files, model checkpoints, etc.
    mkdir -p "$1/results"        # Output files, model checkpoints, etc.

    # Environment file for reproducibility (e.g., Conda)
    touch "$1/environment.yml"

    # Readme file with an overview of the analysis
    touch "$1/README.md"
    echo "Analysis project structure created in '$1'."
}

# Check if an argument (project name) is provided
if [ -z "$1" ]; then
    # If no argument is provided, display an error message
    echo "Error: Please provide a project name."
    echo "Usage: $0 <project_name>"
    exit 1
fi

# Create the analysis folder structure
create_analysis_folder "$1"
