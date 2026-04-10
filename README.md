# af3-ppi
Pipeline for detecting protein-protein interfaces using AlfaFold3

## Installation

### Option 1: Conda (recommended for bioinformatics)

Create and activate a conda environment:

```bash
conda env create -f environment.yml
conda activate af3ppi
```

Install the local `af3ppi` package into the activated environment:

```bash
pip install -e .
```

### Option 2: Virtual environment + pip

Create a virtual environment and install dependencies:

```bash
python3 -m venv af3ppi_env
source af3ppi_env/bin/activate  # On Windows: af3ppi_env\Scripts\activate
python3 -m pip install -r requirements.txt
```

### Option 3: Pip only (if you must)

If you have pip available and want to install globally:

```bash
python3 -m pip install --user -r requirements.txt
```

### Option 4: Install from GitHub

For the latest version directly from GitHub:

```bash
# Install directly from GitHub (creates af3ppi command globally)
pip install git+https://github.com/armstrong-lab/af3-ppi.git

# Or clone and install locally
git clone https://github.com/armstrong-lab/af3-ppi.git
cd af3-ppi
pip install -e .  # Installs af3ppi command
```


Or if you have permission to modify the system Python:

```bash
python3 -m pip install -r requirements.txt
```

## af3-ppi CLI usage

After installation, the package provides a `af3ppi` command. All commands below assume you're in a directory with access to your config files and data.

Generate AF3 JSON AF3_server_inputs:

```bash
af3ppi generate --config config/run_config.yaml --mode binary
```

Binary mode generates pairwise comparisons between every bait and every target in the config file, creating separate AF3 jobs for each bait-target pair.

```bash
af3ppi generate --config config/run_config.yaml --mode multi --target-name PAF1c
```

Multi mode compares each bait to all targets simultaneously in separate analyses, with each job containing one bait and all targets.

```bash
af3ppi generate --config config/run_config.yaml --mode complex
```

Complex mode runs one analysis with all baits and targets combined into a single AF3 job.

## Sequence input types

The config file supports five sequence input styles:

- `whole_proteins`: specify full-length protein names and load the complete protein sequences from the FASTA database.
- `segments`: define protein segments with `protein`, `name`, `start`, and `end`. Only the requested subsequence is used.
- `fasta`: provide custom sequences directly in the YAML file using `name` and `sequence` fields.
- `txt_file`: load a list of sequence names from a text file, where each line contains a protein name.
- `overlapping_windows`: split a protein sequence into multiple overlapping windows. This is useful for scanning large proteins or disordered regions in smaller chunks.

Parse AF3 server outputs and write a summary TSV:

```bash
af3ppi parse-outputs --folder path/to/AF3_server_RNAPII_outputs --out-file results/AF3_server_RNAPII_outputs.txt
```

Generate a single heatmap image from binary output JSONs:

```bash
af3ppi heatmap --config config/run_config.yaml --folder path/to/AF3_server_RNAPII_outputs --out-file results/heatmap_RNAPII.png
```

Generate multi-chain heatmaps from multi output JSONs:

```bash
af3ppi heatmap-multi --config config/run_config.yaml --folder path/to/AF3_server_RNAPII_outputs --out-file results/heatmap_multi_RNAPII.png
```

*Creates a single image file with all heatmaps arranged vertically (one on top of the other)*


## Input types

On the YAML config file, there are five ways to specify the amino acid sequence in the bait or target section:

1 - "whole_protein" - Official gene name. The sequence will be parsed from the FASTA file specified on the "fasta_database_file" field. A default Uniprot library in included.

2 - "segments" - A segment from between the specified start and end of protein defined by the official gene name will be retrived.

3 - "fasta" - Direct input of the sequence

4 - "txt_file" - A text file with the official gene names, one per line, will be read.

5 - "overlapping_windows" - User define overlapping segments of specified start, end, size and overlap.