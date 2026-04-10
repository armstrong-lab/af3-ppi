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
pip install git+https://github.com/armstrong-lab/af3ppi.git

# Or clone and install locally
git clone https://github.com/armstrong-lab/af3ppi.git
cd af3ppi
pip install -e .  # Installs af3ppi command
```


Or if you have permission to modify the system Python:

```bash
python3 -m pip install -r requirements.txt
```

## AF3 CLI usage

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


