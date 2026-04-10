import argparse
from pathlib import Path

from .af3_input import (
    generate_af3_inputs,
    make_output_heatmap,
    make_output_heatmap_multi,
    read_AF3_server_outputs,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="AF3 package CLI for AF3_server_inputs generation, output parsing, and heatmap creation."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    generate_parser = subparsers.add_parser("generate", help="Generate AF3 server JSON AF3_server_inputs.")
    generate_parser.add_argument(
        "--config",
        required=True,
        help="Path to the run configuration YAML file, relative to current directory or absolute.",
    )
    generate_parser.add_argument(
        "--mode",
        choices=["binary", "multi", "complex"],
        required=True,
        help="Type of AF3_server_inputs to generate.",
    )
    generate_parser.add_argument(
        "--target-name",
        help="Target name used for multi mode naming and grouping.",
    )
    generate_parser.add_argument(
        "--output",
        help="Optional output path for the generated JSON AF3_server_inputs.",
    )
    generate_parser.add_argument(
        "--repo-root",
        help="Optional repository root directory. Defaults to the current working directory.",
    )

    parse_parser = subparsers.add_parser(
        "parse-outputs", help="Parse AF3 server output JSONs and write a summary TSV."
    )
    parse_parser.add_argument(
        "--folder",
        required=True,
        help="Path to the folder containing AF3 output JSONs.",
    )
    parse_parser.add_argument(
        "--out-file",
        required=True,
        help="Path to write the summary TSV file.",
    )
    parse_parser.add_argument(
        "--repo-root",
        help="Optional repository root directory. Defaults to the current working directory.",
    )

    heatmap_parser = subparsers.add_parser(
        "heatmap", help="Generate a heatmap image from AF3 summary output JSONs."
    )
    heatmap_parser.add_argument(
        "--config",
        required=True,
        help="Path to the run configuration YAML file used to define bait/target names.",
    )
    heatmap_parser.add_argument(
        "--folder",
        required=True,
        help="Path to the folder containing AF3 output JSONs.",
    )
    heatmap_parser.add_argument(
        "--out-file",
        help="Optional path to save the heatmap image.",
    )
    heatmap_parser.add_argument(
        "--repo-root",
        help="Optional repository root directory. Defaults to the current working directory.",
    )
    heatmap_parser.add_argument(
        "--show",
        action="store_true",
        help="Display the heatmap after generating it.",
    )

    heatmap_multi_parser = subparsers.add_parser(
        "heatmap-multi",
        help="Generate heatmap images for multi-chain AF3 outputs.",
    )
    heatmap_multi_parser.add_argument(
        "--config",
        required=True,
        help="Path to the run configuration YAML file used to define bait/target names.",
    )
    heatmap_multi_parser.add_argument(
        "--folder",
        required=True,
        help="Path to the folder containing AF3 output JSONs.",
    )
    heatmap_multi_parser.add_argument(
        "--out-file",
        help="Optional path to save the combined multi heatmap image.",
    )
    heatmap_multi_parser.add_argument(
        "--repo-root",
        help="Optional repository root directory. Defaults to the current working directory.",
    )
    heatmap_multi_parser.add_argument(
        "--show",
        action="store_true",
        help="Display each heatmap after generating it.",
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    root = Path(args.repo_root) if args.repo_root else Path.cwd()

    if args.command == "generate":
        output_path = generate_af3_inputs(
            config_path=args.config,
            mode=args.mode,
            target_name=args.target_name,
            output_path=args.output,
            root=root,
        )
        print(f"Wrote AF3_server_inputs: {Path(output_path).resolve()}")
        return 0

    if args.command == "parse-outputs":
        output_path = read_AF3_server_outputs(
            folder_name=args.folder,
            out_fname=args.out_file,
            root=root,
        )
        print(f"Wrote summary TSV: {Path(output_path).resolve()}")
        return 0

    if args.command == "heatmap":
        from .af3_input import parse_config, parse_fasta

        run_config = parse_config(args.config, root=root)
        prot_dict = parse_fasta(run_config, root=root)
        output_path = make_output_heatmap(
            prot_dict=prot_dict,
            run_config=run_config,
            output_folder_name=args.folder,
            root=root,
            save_path=args.out_file,
            show=args.show,
        )
        print(f"Wrote heatmap: {Path(output_path).resolve()}")
        return 0

    if args.command == "heatmap-multi":
        from .af3_input import parse_config, parse_fasta

        run_config = parse_config(args.config, root=root)
        prot_dict = parse_fasta(run_config, root=root)
        output_path = make_output_heatmap_multi(
            prot_dict=prot_dict,
            run_config=run_config,
            output_folder_name=args.folder,
            root=root,
            save_path=getattr(args, "out_file", None),
            show=args.show,
        )
        print(f"Wrote combined heatmap: {Path(output_path).resolve()}")
        return 0

    parser.error("Invalid command")


if __name__ == "__main__":
    raise SystemExit(main())
