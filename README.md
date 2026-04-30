# af3-ppi
Pipeline for detecting protein-protein interfaces using AlfaFold3

## Installation

### Download from GitHub

Download the repository from GitHub.



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
- `segments`: define protein segments with `protein`, `name`, `start`, and `end`. Only the requested subsequence is used. You can also use the literal values `start` or `end` to mean the first or last amino acid of the full protein sequence instead of looking up the exact residue number.
- `fasta`: provide custom sequences directly in the YAML file using `name` and `sequence` fields.
- `txt_file`: load a list of sequence names from a text file, where each line contains a protein name.
- `overlapping_windows`: split a protein sequence into multiple overlapping windows. This is useful for scanning large proteins or disordered regions in smaller chunks. For overlapping windows you may also use `start` and `end` to indicate the protein termini when defining the window boundaries.


## Parse AF3 server outputs and write a summary TSV:

```bash
af3ppi parse-outputs --folder path/to/AF3_server_outputs --out-file results/AF3_server_binary_results.txt
```

The summary TSV includes `ipTM`, `minPAE`, and `actifpTM`. The `actifpTM` value is calculated locally from the matching AF3 full data JSON using `pae`, `contact_probs`, and `token_chain_ids`. Interface residues are selected from chain-pair contacts with `contact_probs >= 0.7` by default. Change this with `--contact-cutoff`. The TSV also reports contact support and interface mean pLDDT summaries so small or low-confidence interfaces can be flagged.

```bash
af3ppi parse-outputs --folder path/to/AF3_server_outputs --out-file results/AF3_server_binary_results.txt --contact-cutoff 0.7
```

For multi-chain outputs, write one TSV with one row per bait and separate metric columns for each target:

```bash
af3ppi parse-outputs-multi --config config/run_config.yaml --folder path/to/AF3_server_outputs --out-file results/AF3_server_multi_results.txt
```


Generate a single heatmap image from binary output JSONs:

```bash
af3ppi heatmap --config config/run_config.yaml --folder path/to/AF3_server_outputs --out-file results/heatmap_binary.png
```

Use `--metric` to choose the plotted value:

```bash
af3ppi heatmap --config config/run_config.yaml --folder path/to/AF3_server_outputs --metric actifptm --out-file results/heatmap_binary_actifptm.png
```


Generate multi-chain heatmaps from multi output JSONs:

```bash
af3ppi heatmap-multi --config config/run_config.yaml --folder path/to/AF3_server_outputs --out-file results/heatmap_multi.png
```

The `heatmap-multi` command also accepts `--metric iptm`, `--metric min_pae`, or `--metric actifptm`.

*Creates a single image file with all heatmaps arranged vertically (one on top of the other)*
