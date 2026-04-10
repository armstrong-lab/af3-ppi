import json
import re
from pathlib import Path
from typing import Any, Dict, List

try:
    import yaml
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "PyYAML is required to parse config files. Install it with 'pip install PyYAML'"
    ) from exc

try:
    from Bio import SeqIO
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "Biopython is required to parse FASTA files. Install it with 'pip install biopython'"
    ) from exc


def project_root(root: str | Path | None = None) -> Path:
    if root is not None:
        return Path(root)
    return Path(__file__).resolve().parents[1]


def ensure_parent_dir(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def parse_config(config_fname: str | Path, root: str | Path | None = None) -> Dict[str, Any]:
    repo_root = project_root(root)
    config_path = Path(config_fname)
    if not config_path.is_absolute():
        # If the path already contains 'config', treat it as relative to repo root
        # Otherwise, assume it's a config filename and add 'config/' prefix
        if 'config' in str(config_path):
            config_path = repo_root / config_path
        else:
            config_path = repo_root / "config" / config_path
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    if not config_path.is_file():
        raise ValueError(f"Config path is not a file: {config_path}")
    try:
        with config_path.open("r") as file:
            run_config = yaml.safe_load(file)
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML in {config_path}:\n{exc}") from exc
    if run_config is None:
        raise ValueError(f"Config file is empty: {config_path}")
    
    # Add default run_name if missing
    if 'run_name' not in run_config:
        run_config['run_name'] = config_path.stem
    
    return run_config


def parse_fasta(run_config: Dict[str, Any], root: str | Path | None = None) -> Dict[str, str]:
    repo_root = project_root(root)
    prot_dict: Dict[str, str] = {}
    prot_fname = run_config.get("fasta_database_file")
    if not prot_fname:
        raise KeyError("Missing 'fasta_database_file' key in run_config")
    fasta_path = Path(prot_fname)
    if not fasta_path.is_absolute():
        fasta_path = repo_root / "databases" / fasta_path
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA database file not found: {fasta_path}")
    if not fasta_path.is_file():
        raise ValueError(f"FASTA database path is not a file: {fasta_path}")
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        re_obj = re.search(r"GN=(\S+)", record.description)
        if re_obj:
            prot_dict[re_obj.group(1)] = str(record.seq)
        else:
            prot_dict[record.id] = str(record.seq)
    if not prot_dict:
        raise ValueError(f"No sequences loaded from FASTA file: {fasta_path}")
    return prot_dict


def _load_sequence_by_name(name: str, prot_dict: Dict[str, str]) -> str:
    if name not in prot_dict:
        raise KeyError(f"Protein '{name}' not found in loaded FASTA sequences")
    return prot_dict[name]


def _extract_range(sequence: str, start: str | int, end: str | int) -> str:
    if isinstance(start, str) and start == "start":
        lstart = 0
    else:
        lstart = int(start) - 1
    if isinstance(end, str) and end == "end":
        lend = len(sequence)
    else:
        lend = int(end) - 1
    return sequence[lstart:lend]


def _read_txt_file_names(txt_file: str, root: str | Path | None = None) -> List[str]:
    repo_root = project_root(root)
    txt_path = Path(txt_file)
    if not txt_path.is_absolute():
        txt_path = repo_root / "inputs" / txt_path
    if not txt_path.exists():
        raise FileNotFoundError(f"TXT file not found: {txt_path}")
    if not txt_path.is_file():
        raise ValueError(f"TXT path is not a file: {txt_path}")
    with txt_path.open("r") as file:
        return [line.strip() for line in file if line.strip()]


def get_config_sequences(
    prot_dict: Dict[str, str], key: str, run_config: Dict[str, Any], root: str | Path | None = None
) -> List[Dict[str, str]]:
    if key not in run_config:
        raise KeyError(f"Missing '{key}' section in run_config")
    config_section = run_config[key]
    local_arr: List[Dict[str, str]] = []

    if "segments" in config_section:
        for seg in config_section["segments"]:
            name = seg.get("name")
            if not name:
                raise ValueError("Segment entry missing 'name'")
            sequence = _load_sequence_by_name(seg["protein"], prot_dict)
            seg_seq = _extract_range(sequence, seg["start"], seg["end"])
            local_arr.append({"name": name, "sequence": seg_seq})

    if "whole_proteins" in config_section:
        for protein in config_section["whole_proteins"]:
            name = protein.get("name")
            if not name:
                raise ValueError("Whole protein entry missing 'name'")
            whole_seq = _load_sequence_by_name(name, prot_dict)
            local_arr.append({"name": name, "sequence": whole_seq})

    if "txt_file" in config_section:
        for txt in _read_txt_file_names(config_section["txt_file"], root=root):
            local_arr.append({"name": txt, "sequence": _load_sequence_by_name(txt, prot_dict)})

    if "fasta" in config_section:
        for fasta in config_section["fasta"]:
            if "name" not in fasta or "sequence" not in fasta:
                raise ValueError("FASTA entry must include 'name' and 'sequence'")
            local_arr.append({"name": fasta["name"], "sequence": fasta["sequence"]})

    if "overlapping_windows" in config_section:
        for window in config_section["overlapping_windows"]:
            name = window.get("name")
            if not name:
                raise ValueError("Overlapping window entry missing 'name'")
            main_seq = _load_sequence_by_name(window["protein"], prot_dict)
            start = int(window["start"])
            end = int(window["end"])
            size = int(window["size"])
            overlap = int(window["overlap"])
            if size <= 0 or overlap < 0:
                raise ValueError("Invalid window size/overlap in overlapping_windows")
            if end < start:
                raise ValueError("overlapping_windows end must be >= start")
            lstart = start - 1
            lend = lstart + size
            while lend <= end:
                local_arr.append({"name": f"{name}_{lstart+1}_{lend}", "sequence": main_seq[lstart:lend]})
                lstart = lstart + size - overlap
                lend = lstart + size
            if lstart < end and len(main_seq[lstart:end]) > 0:
                local_arr.append({"name": f"{name}_{lstart+1}_{end}", "sequence": main_seq[lstart:end]})

    if not local_arr:
        raise ValueError(f"No sequences generated for section '{key}' in config")
    return local_arr


def initialize_AF3_server_job(name: str) -> Dict[str, Any]:
    return {
        "name": name,
        "sequences": [],
        "modelSeeds": [],
        "dialect": "alphafoldserver",
        "version": 1,
    }


def add_AF3_job_seq(input_json: Dict[str, Any], seq: str) -> Dict[str, Any]:
    input_json["sequences"].append({"proteinChain": {"sequence": seq, "count": 1}})
    return input_json


def write_manifest(job_arr: List[Dict[str, Any]], output_path: str | Path) -> Path:
    output_path = Path(output_path)
    ensure_parent_dir(output_path)
    with output_path.open("w") as f:
        json.dump(job_arr, f, indent=4)
    return output_path


def make_binary_json_inputs(
    prot_dict: Dict[str, str], run_config: Dict[str, Any], output_path: str | Path | None = None, root: str | Path | None = None
) -> Path:
    bait_arr = get_config_sequences(prot_dict, "bait", run_config, root=root)
    target_arr = get_config_sequences(prot_dict, "target", run_config, root=root)
    if not bait_arr or not target_arr:
        raise ValueError("Bait or target sequences are missing for binary input creation")
    job_arr: List[Dict[str, Any]] = []
    for count, bait in enumerate(bait_arr):
        for target in target_arr:
            tname = f"{count}_{bait['name']}_vs_{target['name']}"
            af3_job = initialize_AF3_server_job(tname)
            add_AF3_job_seq(af3_job, bait["sequence"])
            add_AF3_job_seq(af3_job, target["sequence"])
            job_arr.append(af3_job)
    if output_path is None:
        output_path = project_root(root) / "AF3_server_inputs" / f"AF3_server_binary_{run_config['run_name']}.json"
    return write_manifest(job_arr, output_path)


def make_multi_json_inputs(
    prot_dict: Dict[str, str],
    run_config: Dict[str, Any],
    target_name: str,
    output_path: str | Path | None = None,
    root: str | Path | None = None,
) -> Path:
    bait_arr = get_config_sequences(prot_dict, "bait", run_config, root=root)
    target_arr = get_config_sequences(prot_dict, "target", run_config, root=root)
    if not bait_arr or not target_arr:
        raise ValueError("Bait or target sequences are missing for multi input creation")
    job_arr: List[Dict[str, Any]] = []
    for count, bait in enumerate(bait_arr):
        name = f"{count}_{bait['name']}_{target_name}"
        token_size = len(bait["sequence"])
        af3_job = initialize_AF3_server_job(name)
        add_AF3_job_seq(af3_job, bait["sequence"])
        for target in target_arr:
            add_AF3_job_seq(af3_job, target["sequence"])
            token_size += len(target["sequence"])
        job_arr.append(af3_job)
        print(f"Input: {name}   token size: {token_size}")
    if output_path is None:
        output_path = project_root(root) / "AF3_server_inputs" / f"AF3_server_multi_{run_config['run_name']}.json"
    return write_manifest(job_arr, output_path)


def make_complex_no_bait_input(
    prot_dict: Dict[str, str], run_config: Dict[str, Any], output_path: str | Path | None = None, root: str | Path | None = None
) -> Path:
    lname = run_config.get("run_name")
    if not lname:
        raise KeyError("Missing 'run_name' in run_config")
    af3_job = initialize_AF3_server_job(lname)
    token_size = 0
    if "bait" in run_config:
        for bait in get_config_sequences(prot_dict, "bait", run_config, root=root):
            add_AF3_job_seq(af3_job, bait["sequence"])
            token_size += len(bait["sequence"])
            print(f"Adding {bait['name']}   size: {len(bait['sequence'])}")
    if "target" in run_config:
        for target in get_config_sequences(prot_dict, "target", run_config, root=root):
            add_AF3_job_seq(af3_job, target["sequence"])
            token_size += len(target["sequence"])
            print(f"Adding {target['name']}   size: {len(target['sequence'])}")
    if not af3_job["sequences"]:
        raise ValueError("No sequences added to complex input; check bait/target configuration")
    if output_path is None:
        output_path = project_root(root) / "AF3_server_inputs" / f"{lname}_inputs.json"
    output_path = Path(output_path)
    return write_manifest([af3_job], output_path)


def _ensure_visualization_deps():
    try:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        import numpy as np
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "pandas, seaborn, matplotlib, and numpy are required for output parsing and heatmaps. "
            "Install them with 'pip install pandas seaborn matplotlib numpy'"
        ) from exc
    return pd, sns, plt, np


def read_AF3_server_outputs(folder_name: str, out_fname: str, root: str | Path | None = None) -> Path:
    output_folder = Path(folder_name)
    if not output_folder.is_absolute():
        output_folder = project_root(root) / output_folder

    if not output_folder.exists():
        raise FileNotFoundError(f"AF3 server output folder not found: {output_folder}")
    pd, _, _, _ = _ensure_visualization_deps()
    outdata = pd.DataFrame(columns=["ipTM", "minPEA"])
    found = False
    for f in output_folder.rglob("*summary_confidences_0.json"):
        found = True
        with f.open("r") as file:
            jdata = json.load(file)
        re_obj = re.search(r"fold_(.+)_summary", f.stem)
        if not re_obj:
            continue
        run_name = re_obj.group(1)
        min_pae = min(jdata["chain_pair_pae_min"][0][1], jdata["chain_pair_pae_min"][1][0])
        iptm = jdata["iptm"]
        outdata.loc[run_name] = [iptm, min_pae]
    if not found:
        raise ValueError(f"No summary confidence JSON files found in {output_folder}")
    
    out_path = Path(out_fname)
    if not out_path.is_absolute():
        output_path = project_root(root) / out_path
    else:
        output_path = out_path
    
    ensure_parent_dir(output_path)
    outdata.to_csv(output_path, sep="\t")
    return output_path


# compatibility alias for existing notebook usage
read_AF3sever_outputs = read_AF3_server_outputs


def make_output_heatmap(
    prot_dict: Dict[str, str],
    run_config: Dict[str, Any],
    output_folder_name: str,
    root: str | Path | None = None,
    save_path: str | Path | None = None,
    show: bool = False,
) -> Path:
    pd, sns, plt, np = _ensure_visualization_deps()
    bait_arr = get_config_sequences(prot_dict, "bait", run_config, root=root)
    target_arr = get_config_sequences(prot_dict, "target", run_config, root=root)
    bait_names = [bait["name"].upper() for bait in bait_arr]
    target_names = [target["name"].upper() for target in target_arr]
    df = pd.DataFrame(index=bait_names, columns=target_names, dtype=float)
    output_folder = Path(output_folder_name)
    if not output_folder.is_absolute():
        output_folder = project_root(root) / output_folder
    if not output_folder.exists():
        raise FileNotFoundError(f"AF3 server output folder not found: {output_folder}")
    found = False
    for f in output_folder.rglob("*summary_confidences_0.json"):
        with f.open("r") as file:
            jdata = json.load(file)
        re_obj = re.search(r"fold_\d+_(\S+)_vs_(\S+)_summary", f.stem)
        if not re_obj:
            continue
        found = True
        bait_name = re_obj.group(1).upper()
        target_name = re_obj.group(2).upper()
        df.loc[bait_name, target_name] = float(jdata["iptm"])
    if not found:
        raise ValueError(f"No heatmap source JSON files found in {output_folder}")
    df = df.apply(pd.to_numeric, errors="coerce")
    plt.figure(figsize=(max(8, len(target_names) * 0.5), max(6, len(bait_names) * 0.5)))
    sns.heatmap(df, vmin=0.2, vmax=1)
    if save_path is None:
        # Use folder name for filename to avoid path separators
        folder_name = Path(output_folder_name).name
        save_path = project_root(root) / "results" / f"heatmap_{folder_name}.png"
    save_path = Path(save_path)
    ensure_parent_dir(save_path)
    plt.savefig(save_path, bbox_inches="tight")
    if show:
        plt.show()
    plt.close()
    return save_path


def make_output_heatmap_multi(
    prot_dict: Dict[str, str],
    run_config: Dict[str, Any],
    output_folder_name: str,
    root: str | Path | None = None,
    save_path: str | Path | None = None,
    show: bool = False,
) -> Path:
    pd, sns, plt, np = _ensure_visualization_deps()
    bait_arr = get_config_sequences(prot_dict, "bait", run_config, root=root)
    target_arr = get_config_sequences(prot_dict, "target", run_config, root=root)
    target_name_arr = [target["name"] for target in target_arr]
    output_folder = Path(output_folder_name)
    if not output_folder.is_absolute():
        output_folder = project_root(root) / output_folder
    if not output_folder.exists():
        raise FileNotFoundError(f"AF3 server output folder not found: {output_folder}")
    df_dict: Dict[str, pd.DataFrame] = {}
    found_bait_names: set[str] = set()
    found = False
    for f in output_folder.rglob("*summary_confidences_0.json"):
        re_obj = re.search(r"fold_\d+_(.+)_[^_]+_summary", f.stem)
        if not re_obj:
            continue
        found = True
        bait_name = re_obj.group(1).upper()
        found_bait_names.add(bait_name)
        mat_arr = [bait_name] + target_name_arr
        with f.open("r") as file:
            jdata = json.load(file)
        lmat = np.array(jdata["chain_pair_iptm"])
        if lmat.shape != (len(mat_arr), len(mat_arr)):
            raise ValueError(
                f"Chain pair matrix shape mismatch for file {f}. Expected ({len(mat_arr)}, {len(mat_arr)}) based on config (bait + {len(target_arr)} targets), got {lmat.shape}. This indicates the config does not match the AF3 results."
            )
        df_dict[bait_name] = pd.DataFrame(lmat, index=mat_arr, columns=mat_arr)
    if not found:
        raise ValueError(f"No multi heatmap JSON files found in {output_folder}")
    
    # Validate that config baits match found baits
    config_bait_names = [bait["name"].upper() for bait in bait_arr]
    if not any(cb in found_bait_names for cb in config_bait_names):
        raise ValueError(
            f"Config baits {config_bait_names} do not match any baits found in results {list(found_bait_names)}. Please ensure the config matches the results folder."
        )
    
    # Filter baits that have data
    available_baits = [bait for bait in bait_arr if bait["name"].upper() in df_dict]
    
    if not available_baits:
        raise ValueError("No bait data found for heatmap generation")
    
    # Create single figure with subplots arranged vertically
    n_baits = len(available_baits)
    fig, axes = plt.subplots(n_baits, 1, figsize=(max(8, len(mat_arr) * 0.5), max(6, len(mat_arr) * 0.5) * n_baits))
    
    # Handle single subplot case
    if n_baits == 1:
        axes = [axes]
    
    saved_paths: List[Path] = []
    for i, bait in enumerate(available_baits):
        bait_name = bait["name"].upper()
        sns.heatmap(df_dict[bait_name], vmin=0.2, vmax=1, ax=axes[i])
        axes[i].set_title(bait["name"], fontsize=16, color="blue", loc="center")
    
    plt.tight_layout()
    
    # Save as single combined heatmap file
    if save_path is None:
        # Use the folder name (not full path) for the default save directory
        folder_name = Path(output_folder_name).name if Path(output_folder_name).is_absolute() else Path(output_folder_name).name
        save_dir = project_root(root) / "results"
        output_path = save_dir / f"heatmap_multi_{folder_name}.png"
    else:
        output_path = Path(save_path)
        ensure_parent_dir(output_path)
    
    ensure_parent_dir(output_path)
    plt.savefig(output_path, bbox_inches="tight")
    if show:
        plt.show()
    plt.close()
    
    return output_path


def generate_af3_inputs(
    config_path: str | Path,
    mode: str,
    target_name: str | None = None,
    output_path: str | Path | None = None,
    root: str | Path | None = None,
) -> Path:
    run_config = parse_config(config_path, root=root)
    prot_dict = parse_fasta(run_config, root=root)
    if mode == "binary":
        return make_binary_json_inputs(prot_dict, run_config, output_path=output_path, root=root)
    if mode == "multi":
        if not target_name:
            raise ValueError("target_name is required for multi mode")
        return make_multi_json_inputs(prot_dict, run_config, target_name, output_path=output_path, root=root)
    if mode == "complex":
        return make_complex_no_bait_input(prot_dict, run_config, output_path=output_path, root=root)
    raise ValueError(f"Unsupported mode: {mode}")
