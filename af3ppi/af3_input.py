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
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
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


def _ensure_pandas_dep():
    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "pandas is required for output parsing. Install it with 'pip install pandas'"
        ) from exc
    return pd


def _ensure_numpy_dep():
    try:
        import numpy as np
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "numpy is required for local actifpTM calculation. Install it with 'pip install numpy'"
        ) from exc
    return np


METRIC_ALIASES = {
    "iptm": "iptm",
    "ipTM": "iptm",
    "min_pae": "min_pae",
    "min-pae": "min_pae",
    "minpae": "min_pae",
    "actifptm": "actifptm",
    "actifpTM": "actifptm",
}


def normalize_metric(metric: str) -> str:
    metric_key = METRIC_ALIASES.get(metric, METRIC_ALIASES.get(metric.lower()))
    if metric_key is None:
        valid_metrics = ", ".join(sorted({"iptm", "min_pae", "actifptm"}))
        raise ValueError(f"Unsupported metric '{metric}'. Choose one of: {valid_metrics}")
    return metric_key


def metric_display_name(metric: str) -> str:
    metric = normalize_metric(metric)
    return {"iptm": "ipTM", "min_pae": "minPAE", "actifptm": "actifpTM"}[metric]


def _get_symmetric_pair(matrix: Any, i: int, j: int) -> float:
    return float(max(matrix[i][j], matrix[j][i]))


ACTIFPTM_CONTACT_PROB_CUTOFF = 0.7


def validate_contact_prob_cutoff(contact_prob_cutoff: float) -> float:
    contact_prob_cutoff = float(contact_prob_cutoff)
    if contact_prob_cutoff < 0 or contact_prob_cutoff > 1:
        raise ValueError("actifpTM contact probability cutoff must be between 0 and 1")
    return contact_prob_cutoff


def _get_pairwise_dict_value(pairwise: Dict[str, Any], i: int, j: int) -> float:
    chain_i = chr(ord("A") + i)
    chain_j = chr(ord("A") + j)
    for key in (f"{chain_i}-{chain_j}", f"{chain_j}-{chain_i}", f"{i}-{j}", f"{j}-{i}"):
        if key in pairwise:
            return float(pairwise[key])
    raise KeyError(f"No pairwise metric entry found for chains {chain_i}-{chain_j}")


def _find_full_confidence_json(summary_path: Path) -> Path:
    candidates = [
        summary_path.with_name(summary_path.name.replace("summary_confidences", "full_data")),
        summary_path.with_name(summary_path.name.replace("summary_confidences", "confidences")),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Could not find full AF3 confidence JSON for {summary_path}. Expected a matching "
        "'full_data' or 'confidences' JSON in the same folder."
    )


def _find_model_cif(summary_path: Path) -> Path | None:
    candidates = [
        summary_path.with_name(summary_path.name.replace("summary_confidences", "model")).with_suffix(".cif"),
        summary_path.with_name(re.sub(r"_summary_confidences(?:_\d+)?\.json$", "_model.cif", summary_path.name)),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _load_full_confidence_json(summary_path: Path | None) -> Dict[str, Any]:
    if summary_path is None:
        return {}
    full_path = _find_full_confidence_json(summary_path)
    with full_path.open("r") as file:
        return json.load(file)


def _chain_token_indices(full_data: Dict[str, Any]) -> List[Any]:
    if "token_chain_ids" not in full_data:
        raise KeyError("Missing token_chain_ids in full confidence JSON")
    np = _ensure_numpy_dep()
    chain_ids = np.array(full_data["token_chain_ids"])
    return [np.where(chain_ids == chain_id)[0] for chain_id in dict.fromkeys(chain_ids.tolist())]


def _as_list(value: Any) -> List[Any]:
    return value if isinstance(value, list) else [value]


def _token_plddt_map_from_cif(cif_path: Path | None) -> Dict[tuple[str, int], float]:
    if cif_path is None:
        return {}
    cif_data = MMCIF2Dict(str(cif_path))
    chain_ids = _as_list(
        cif_data.get("_atom_site.auth_asym_id")
        or cif_data.get("_atom_site.label_asym_id")
        or []
    )
    res_ids = _as_list(
        cif_data.get("_atom_site.auth_seq_id")
        or cif_data.get("_atom_site.label_seq_id")
        or []
    )
    b_factors = _as_list(cif_data.get("_atom_site.B_iso_or_equiv") or [])
    if not chain_ids or not res_ids or not b_factors:
        return {}

    residue_values: Dict[tuple[str, int], List[float]] = {}
    for chain_id, res_id, b_factor in zip(chain_ids, res_ids, b_factors):
        if res_id in {".", "?"}:
            continue
        try:
            key = (str(chain_id), int(float(res_id)))
            residue_values.setdefault(key, []).append(float(b_factor))
        except ValueError:
            continue
    return {
        key: sum(values) / len(values)
        for key, values in residue_values.items()
        if values
    }


def _token_plddts(full_data: Dict[str, Any], token_indices: Any, plddt_map: Dict[tuple[str, int], float]) -> List[float]:
    if "token_plddts" in full_data:
        return [float(full_data["token_plddts"][int(idx)]) for idx in token_indices]
    if len(full_data.get("atom_plddts", [])) == len(full_data.get("token_chain_ids", [])):
        return [float(full_data["atom_plddts"][int(idx)]) for idx in token_indices]
    if not plddt_map:
        return []
    token_chain_ids = full_data.get("token_chain_ids", [])
    token_res_ids = full_data.get("token_res_ids", [])
    plddts = []
    for idx in token_indices:
        token_idx = int(idx)
        key = (str(token_chain_ids[token_idx]), int(token_res_ids[token_idx]))
        if key in plddt_map:
            plddts.append(float(plddt_map[key]))
    return plddts


def _mean_or_nan(values: List[float]) -> float:
    if not values:
        return float("nan")
    return round(sum(values) / len(values), 2)


def calculate_actifptm_from_af3(
    full_data: Dict[str, Any],
    i: int = 0,
    j: int = 1,
    contact_prob_cutoff: float = ACTIFPTM_CONTACT_PROB_CUTOFF,
) -> float:
    contact_prob_cutoff = validate_contact_prob_cutoff(contact_prob_cutoff)
    return calculate_actifptm_details_from_af3(
        full_data,
        i=i,
        j=j,
        contact_prob_cutoff=contact_prob_cutoff,
    )["actifptm"]


def calculate_actifptm_details_from_af3(
    full_data: Dict[str, Any],
    i: int = 0,
    j: int = 1,
    contact_prob_cutoff: float = ACTIFPTM_CONTACT_PROB_CUTOFF,
    plddt_map: Dict[tuple[str, int], float] | None = None,
) -> Dict[str, float | int]:
    contact_prob_cutoff = validate_contact_prob_cutoff(contact_prob_cutoff)
    np = _ensure_numpy_dep()
    if "pae" not in full_data:
        raise KeyError("Missing pae in full confidence JSON")
    if "contact_probs" not in full_data:
        raise KeyError("Missing contact_probs in full confidence JSON")

    pae = np.array(full_data["pae"], dtype=float)
    contact_probs = np.array(full_data["contact_probs"], dtype=float)
    chain_indices = _chain_token_indices(full_data)
    if i >= len(chain_indices) or j >= len(chain_indices):
        raise IndexError(f"Chain pair ({i}, {j}) is outside the {len(chain_indices)} chains in token_chain_ids")

    idx_i = chain_indices[i]
    idx_j = chain_indices[j]
    pair_contacts = contact_probs[np.ix_(idx_i, idx_j)] >= contact_prob_cutoff
    contact_count = int(pair_contacts.sum())
    if not pair_contacts.any():
        return {
            "actifptm": 0.0,
            "contact_count": 0,
            "chain_i_interface_residues": 0,
            "chain_j_interface_residues": 0,
            "interface_mean_plddt": float("nan"),
            "chain_i_interface_mean_plddt": float("nan"),
            "chain_j_interface_mean_plddt": float("nan"),
        }

    interface_i = idx_i[np.where(pair_contacts)[0]]
    interface_j = idx_j[np.where(pair_contacts)[1]]
    unique_interface_i = np.unique(interface_i)
    unique_interface_j = np.unique(interface_j)
    chain_i_interface_residues = int(len(np.unique(interface_i)))
    chain_j_interface_residues = int(len(np.unique(interface_j)))
    selected = np.unique(np.concatenate((unique_interface_i, unique_interface_j)))
    if plddt_map is None:
        plddt_map = {}
    chain_i_plddts = _token_plddts(full_data, unique_interface_i, plddt_map)
    chain_j_plddts = _token_plddts(full_data, unique_interface_j, plddt_map)
    interface_plddts = chain_i_plddts + chain_j_plddts
    residue_weights = np.zeros(pae.shape[0], dtype=float)
    residue_weights[selected] = 1.0

    pair_mask = np.zeros(pae.shape, dtype=float)
    pair_mask[np.ix_(idx_i, idx_j)] = 1.0
    pair_mask[np.ix_(idx_j, idx_i)] = 1.0
    pair_residue_weights = pair_mask * residue_weights[:, None] * residue_weights[None, :]

    num_res = pae.shape[0]
    clipped_num_res = max(num_res, 19)
    d0 = 1.24 * (clipped_num_res - 15) ** (1.0 / 3.0) - 1.8
    predicted_tm_term = 1.0 / (1.0 + np.square(pae) / np.square(d0))
    normed_weights = pair_residue_weights / (1e-8 + pair_residue_weights.sum(-1, keepdims=True))
    residuewise_actifptm = (predicted_tm_term * normed_weights).sum(-1) * residue_weights
    return {
        "actifptm": round(float(residuewise_actifptm.max()), 3),
        "contact_count": contact_count,
        "chain_i_interface_residues": chain_i_interface_residues,
        "chain_j_interface_residues": chain_j_interface_residues,
        "interface_mean_plddt": _mean_or_nan(interface_plddts),
        "chain_i_interface_mean_plddt": _mean_or_nan(chain_i_plddts),
        "chain_j_interface_mean_plddt": _mean_or_nan(chain_j_plddts),
    }


def get_actifptm_support_metrics(
    summary_path: Path,
    i: int = 0,
    j: int = 1,
    contact_prob_cutoff: float = ACTIFPTM_CONTACT_PROB_CUTOFF,
) -> Dict[str, float | int]:
    full_data = _load_full_confidence_json(summary_path)
    return calculate_actifptm_details_from_af3(
        full_data,
        i=i,
        j=j,
        contact_prob_cutoff=contact_prob_cutoff,
        plddt_map=_token_plddt_map_from_cif(_find_model_cif(summary_path)),
    )


def calculate_actifptm_matrix_from_af3(
    full_data: Dict[str, Any],
    expected_size: int,
    np: Any,
    summary_data: Dict[str, Any] | None = None,
    contact_prob_cutoff: float = ACTIFPTM_CONTACT_PROB_CUTOFF,
) -> Any:
    matrix = np.zeros((expected_size, expected_size), dtype=float)
    if summary_data is not None and "chain_ptm" in summary_data:
        for i, chain_ptm in enumerate(summary_data["chain_ptm"][:expected_size]):
            matrix[i, i] = float(chain_ptm)
    for i in range(expected_size):
        for j in range(i + 1, expected_size):
            value = calculate_actifptm_from_af3(
                full_data,
                i=i,
                j=j,
                contact_prob_cutoff=contact_prob_cutoff,
            )
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix


def get_metric_value(
    jdata: Dict[str, Any],
    metric: str,
    i: int = 0,
    j: int = 1,
    summary_path: Path | None = None,
    contact_prob_cutoff: float = ACTIFPTM_CONTACT_PROB_CUTOFF,
) -> float:
    metric = normalize_metric(metric)
    if metric == "iptm":
        if "iptm" in jdata:
            return float(jdata["iptm"])
        if "chain_pair_iptm" in jdata:
            return _get_symmetric_pair(jdata["chain_pair_iptm"], i, j)
        if "pairwise_iptm" in jdata and isinstance(jdata["pairwise_iptm"], dict):
            return _get_pairwise_dict_value(jdata["pairwise_iptm"], i, j)
        raise KeyError("Missing iptm, chain_pair_iptm, or pairwise_iptm in summary JSON")
    if metric == "min_pae":
        return _get_symmetric_pair(jdata["chain_pair_pae_min"], i, j)
    if metric == "actifptm":
        for matrix_key in ("chain_pair_actifptm", "chain_pair_actifpTM", "pairwise_actifptm_matrix"):
            if matrix_key in jdata:
                return _get_symmetric_pair(jdata[matrix_key], i, j)
        if "pairwise_actifptm" in jdata and isinstance(jdata["pairwise_actifptm"], dict):
            return _get_pairwise_dict_value(jdata["pairwise_actifptm"], i, j)
        for scalar_key in ("actifptm", "actifpTM"):
            if scalar_key in jdata:
                return float(jdata[scalar_key])
        return calculate_actifptm_from_af3(
            _load_full_confidence_json(summary_path),
            i=i,
            j=j,
            contact_prob_cutoff=contact_prob_cutoff,
        )
    raise ValueError(f"Unsupported metric '{metric}'")


def get_metric_matrix(
    jdata: Dict[str, Any],
    metric: str,
    expected_size: int,
    np: Any,
    summary_path: Path | None = None,
    contact_prob_cutoff: float = ACTIFPTM_CONTACT_PROB_CUTOFF,
) -> Any:
    metric = normalize_metric(metric)
    if metric == "iptm":
        if "chain_pair_iptm" in jdata:
            return np.array(jdata["chain_pair_iptm"], dtype=float)
        if "pairwise_iptm" in jdata and isinstance(jdata["pairwise_iptm"], dict):
            return _pairwise_dict_to_matrix(jdata["pairwise_iptm"], expected_size, np)
        raise KeyError("Missing chain_pair_iptm or pairwise_iptm in summary JSON")
    if metric == "min_pae":
        return np.array(jdata["chain_pair_pae_min"], dtype=float)
    if metric == "actifptm":
        for matrix_key in ("chain_pair_actifptm", "chain_pair_actifpTM", "pairwise_actifptm_matrix"):
            if matrix_key in jdata:
                return np.array(jdata[matrix_key], dtype=float)
        if "pairwise_actifptm" in jdata and isinstance(jdata["pairwise_actifptm"], dict):
            matrix = _pairwise_dict_to_matrix(jdata["pairwise_actifptm"], expected_size, np)
            per_chain = jdata.get("per_chain_ptm")
            if isinstance(per_chain, dict):
                for i in range(expected_size):
                    chain_label = chr(ord("A") + i)
                    if chain_label in per_chain:
                        matrix[i, i] = float(per_chain[chain_label])
            return matrix
        return calculate_actifptm_matrix_from_af3(
            _load_full_confidence_json(summary_path),
            expected_size,
            np,
            summary_data=jdata,
            contact_prob_cutoff=contact_prob_cutoff,
        )
    raise ValueError(f"Unsupported metric '{metric}'")


def _pairwise_dict_to_matrix(pairwise: Dict[str, Any], expected_size: int, np: Any) -> Any:
    matrix = np.full((expected_size, expected_size), np.nan, dtype=float)
    for key, value in pairwise.items():
        re_obj = re.fullmatch(r"([A-Z]|\d+)-([A-Z]|\d+)", str(key))
        if not re_obj:
            continue
        left, right = re_obj.groups()
        i = ord(left) - ord("A") if left.isalpha() else int(left)
        j = ord(right) - ord("A") if right.isalpha() else int(right)
        if i >= expected_size or j >= expected_size:
            continue
        matrix[i, j] = float(value)
        matrix[j, i] = float(value)
    return matrix


def _heatmap_kwargs(metric: str) -> Dict[str, Any]:
    metric = normalize_metric(metric)
    if metric == "min_pae":
        return {"cmap": "viridis_r"}
    return {"vmin": 0.2, "vmax": 1}


def read_AF3_server_outputs(
    folder_name: str,
    out_fname: str,
    root: str | Path | None = None,
    actifptm_contact_prob_cutoff: float = ACTIFPTM_CONTACT_PROB_CUTOFF,
) -> Path:
    actifptm_contact_prob_cutoff = validate_contact_prob_cutoff(actifptm_contact_prob_cutoff)
    output_folder = Path(folder_name)
    if not output_folder.is_absolute():
        output_folder = project_root(root) / output_folder

    if not output_folder.exists():
        raise FileNotFoundError(f"AF3 server output folder not found: {output_folder}")
    pd = _ensure_pandas_dep()
    out_rows: List[Dict[str, Any]] = []
    found = False
    for f in output_folder.rglob("*summary_confidences_0.json"):
        found = True
        with f.open("r") as file:
            jdata = json.load(file)
        re_obj = re.search(r"fold_(.+)_summary", f.stem)
        if not re_obj:
            continue
        run_name = re_obj.group(1)
        iptm = get_metric_value(jdata, "iptm")
        min_pae = get_metric_value(jdata, "min_pae")
        actifptm_details = get_actifptm_support_metrics(
            f,
            contact_prob_cutoff=actifptm_contact_prob_cutoff,
        )
        out_rows.append(
            {
                "run_name": run_name,
                "ipTM": iptm,
                "minPAE": min_pae,
                "actifpTM": actifptm_details["actifptm"],
                "actifpTM_contacts": actifptm_details["contact_count"],
                "actifpTM_chain_i_residues": actifptm_details["chain_i_interface_residues"],
                "actifpTM_chain_j_residues": actifptm_details["chain_j_interface_residues"],
                "actifpTM_interface_mean_pLDDT": actifptm_details["interface_mean_plddt"],
                "actifpTM_chain_i_mean_pLDDT": actifptm_details["chain_i_interface_mean_plddt"],
                "actifpTM_chain_j_mean_pLDDT": actifptm_details["chain_j_interface_mean_plddt"],
            }
        )
    if not found:
        raise ValueError(f"No summary confidence JSON files found in {output_folder}")
    if not out_rows:
        raise ValueError(f"No parseable summary confidence JSON files found in {output_folder}")
    
    out_path = Path(out_fname)
    if not out_path.is_absolute():
        output_path = project_root(root) / out_path
    else:
        output_path = out_path
    
    ensure_parent_dir(output_path)
    outdata = pd.DataFrame(out_rows).set_index("run_name")
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
    metric: str = "iptm",
    actifptm_contact_prob_cutoff: float = ACTIFPTM_CONTACT_PROB_CUTOFF,
) -> Path:
    pd, sns, plt, np = _ensure_visualization_deps()
    metric = normalize_metric(metric)
    actifptm_contact_prob_cutoff = validate_contact_prob_cutoff(actifptm_contact_prob_cutoff)
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
        df.loc[bait_name, target_name] = get_metric_value(
            jdata,
            metric,
            summary_path=f,
            contact_prob_cutoff=actifptm_contact_prob_cutoff,
        )
    if not found:
        raise ValueError(f"No heatmap source JSON files found in {output_folder}")
    df = df.apply(pd.to_numeric, errors="coerce")
    plt.figure(figsize=(max(8, len(target_names) * 0.5), max(6, len(bait_names) * 0.5)))
    sns.heatmap(df, **_heatmap_kwargs(metric))
    plt.title(metric_display_name(metric))
    if save_path is None:
        # Use folder name for filename to avoid path separators
        folder_name = Path(output_folder_name).name
        save_path = project_root(root) / "results" / f"heatmap_{metric}_{folder_name}.png"
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
    metric: str = "iptm",
    actifptm_contact_prob_cutoff: float = ACTIFPTM_CONTACT_PROB_CUTOFF,
) -> Path:
    pd, sns, plt, np = _ensure_visualization_deps()
    metric = normalize_metric(metric)
    actifptm_contact_prob_cutoff = validate_contact_prob_cutoff(actifptm_contact_prob_cutoff)
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
        lmat = get_metric_matrix(
            jdata,
            metric,
            len(mat_arr),
            np,
            summary_path=f,
            contact_prob_cutoff=actifptm_contact_prob_cutoff,
        )
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
        sns.heatmap(df_dict[bait_name], ax=axes[i], **_heatmap_kwargs(metric))
        axes[i].set_title(f"{bait['name']} - {metric_display_name(metric)}", fontsize=16, color="blue", loc="center")
    
    plt.tight_layout()
    
    # Save as single combined heatmap file
    if save_path is None:
        # Use the folder name (not full path) for the default save directory
        folder_name = Path(output_folder_name).name if Path(output_folder_name).is_absolute() else Path(output_folder_name).name
        save_dir = project_root(root) / "results"
        output_path = save_dir / f"heatmap_multi_{metric}_{folder_name}.png"
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
