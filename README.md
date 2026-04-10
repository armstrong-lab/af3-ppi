# af3-ppi
Pipeline for detecting protein-protein interfaces using AlfaFold3

## Installation

### Download from GitHub

Download the repository from GitHub:



### Set up environment (recommended)

Create and activate a conda environment:

```bash
conda env create -f environment.yml
conda activate af3ppi
```

### Alternative: Install dependencies with pip

If you prefer not to use conda, install dependencies with pip:

```bash
pip install -r requirements.txt
```

## af3-ppi CLI usage

After setting up the environment, run the tool from the repository directory:

```bash
python -m af3ppi --help
```

All commands below assume you're in the repository directory with the environment activated.

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


