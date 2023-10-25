# Scanpy to CellRanger MTX Format Converter

This script provides a utility to convert a Scanpy annotated data matrix (`adata`) into the CellRanger v3 MTX format. The output will be saved in the specified directory and will contain three files: `features.tsv.gz`, `barcodes.tsv.gz`, and `matrix.mtx.gz`.

## Dependencies:
- pandas
- scipy
- os, gzip (standard libraries)

## Usage:

Firstly, ensure that you have the required libraries installed:

```bash
pip install pandas scipy
```

Next, import the function `write_mtx` from the provided script and use it in your Python environment:

```python
from your_script_name import write_mtx
import scanpy as sc

# Load your data into an AnnData object (this is just a general example)
adata = sc.read_10x_h5('path_to_your_data.h5')

# Use the function to convert and save
output_directory = 'desired_output_directory'
write_mtx(adata, output_directory)
```

Replace `'path_to_your_data.h5'` with the path to your dataset and `'desired_output_directory'` with the directory where you want the output to be saved.

## Output:

The script will create (or overwrite if they exist) the following files in the specified output directory:

- `features.tsv.gz`: Contains the gene identifiers.
- `barcodes.tsv.gz`: Contains the cell barcodes.
- `matrix.mtx.gz`: Contains the gene expression values in the MTX format.

---

To use this README, save the content into a `README.md` file in the same directory as your script. You can then view it with any Markdown viewer or on platforms like GitHub which automatically render Markdown files.
