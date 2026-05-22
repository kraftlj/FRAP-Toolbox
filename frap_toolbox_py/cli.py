from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

from .data.loading import load_diffusion_datasets, load_reaction_datasets
from .exports import write_analysis_bundle
from .models.diffusion import default_diffusion_config, fit_diffusion_model
from .models.reaction import (
    default_reaction1_config,
    default_reaction2_config,
    fit_reaction_model,
)
from .roi import CircularROI, adjacent_circle, load_roi_mask
from .types import BasicInputs


def _build_mask_factory(roi: CircularROI):
    def factory(shape):
        return roi.to_mask(shape)

    return factory


def _load_mask_argument(path: Path, name: str):
    try:
        return load_roi_mask(path, name, allow_single=True)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"Could not load {name} ROI mask from {path}: {exc}") from exc


def parse_arguments(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="FRAP-Toolbox Python port")
    parser.add_argument("files", nargs="+", type=Path, help="Paths to FRAP image files")
    parser.add_argument("--model", choices=["diffusion", "reaction1", "reaction2"], default="diffusion")
    parser.add_argument(
        "--roi",
        nargs=3,
        type=float,
        metavar=("X", "Y", "RADIUS"),
        help="Circular bleach ROI definition in pixels",
    )
    parser.add_argument(
        "--bleach-mask",
        type=Path,
        help="Saved .npz ROI mask file containing a 'bleach' mask",
    )
    parser.add_argument("--post-bleach-frame", type=int, default=21)
    parser.add_argument("--pre-bleach-count", type=int, default=10)
    parser.add_argument("--background", type=float, default=0.0)
    parser.add_argument("--normalize-by-cell", action="store_true")
    parser.add_argument(
        "--cell-roi",
        nargs=3,
        type=float,
        metavar=("X", "Y", "RADIUS"),
        help="Circular whole-cell ROI definition in pixels for reaction whole-cell normalization",
    )
    parser.add_argument(
        "--cell-mask",
        type=Path,
        help="Saved .npz ROI mask file containing a 'cell' mask",
    )
    parser.add_argument("--use-adjacent-roi", action="store_true")
    parser.add_argument(
        "--adjacent-offset",
        type=float,
        default=2.5,
        help="Offset multiplier for adjacent circular ROI (used with --use-adjacent-roi)",
    )
    parser.add_argument(
        "--fit-mode",
        choices=[
            "global",
            "individual",
            "average_curve",
            "simplified_kang",
            "simplified_kang_global",
        ],
        default="global",
        help=(
            "Diffusion fit strategy. 'global' concatenates per-curve residuals without averaging curves; "
            "'simplified_kang' uses the Kang half-time estimator; 'simplified_kang_global' fits the "
            "simplified recovery equation to pooled curves. For reaction models, 'individual' and "
            "'average_curve' select per-curve or averaged-curve fitting."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Write the local beta export bundle to this directory.",
    )
    return parser.parse_args(argv)


def _image_shape_from_datasets(datasets) -> tuple[int, int] | None:
    for dataset in datasets:
        image_shape = dataset.metadata.get("image_shape")
        if image_shape is not None:
            return tuple(int(value) for value in image_shape)
    return None


def _mask_for_export(mask_or_factory, image_shape: tuple[int, int] | None):
    if mask_or_factory is None or image_shape is None:
        return None
    if hasattr(mask_or_factory, "shape"):
        return mask_or_factory
    return mask_or_factory(image_shape)


def _write_cli_exports(args, result, datasets, mask_factory, cell_factory, model: str) -> None:
    if args.output_dir is None:
        return

    image_shape = _image_shape_from_datasets(datasets)
    roi_masks = {}
    bleach_mask = _mask_for_export(mask_factory, image_shape)
    if bleach_mask is not None:
        roi_masks["bleach"] = bleach_mask
    cell_mask = _mask_for_export(cell_factory, image_shape)
    if cell_mask is not None:
        roi_masks["cell"] = cell_mask

    roi_sources = {
        "bleach": "Saved mask upload" if args.bleach_mask is not None else "Circular numeric ROI",
    }
    if args.cell_mask is not None:
        roi_sources["cell"] = "Saved mask upload"
    elif args.cell_roi is not None:
        roi_sources["cell"] = "Circular numeric ROI"

    write_analysis_bundle(
        args.output_dir,
        result,
        model=model,
        files=args.files,
        settings={
            "post_bleach_frame": args.post_bleach_frame,
            "pre_bleach_count": args.pre_bleach_count,
            "background": args.background,
            "normalize_by_cell": args.normalize_by_cell,
            "use_adjacent_roi": args.use_adjacent_roi,
            "adjacent_offset": args.adjacent_offset,
        },
        roi_sources=roi_sources,
        fit_mode=args.fit_mode,
        roi_masks=roi_masks,
        roi_extra_metadata={
            "frame_index": None,
            "roi_source": roi_sources,
            "model": model,
            "created_by": "frap-toolbox-cli",
        },
        created_by="frap-toolbox-cli",
    )


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_arguments(argv)

    if args.roi is None and args.bleach_mask is None:
        raise SystemExit("Provide either --roi X Y RADIUS or --bleach-mask PATH.")
    if args.roi is not None and args.bleach_mask is not None:
        raise SystemExit("Use either --roi or --bleach-mask for the bleach ROI, not both.")
    if args.cell_roi is not None and args.cell_mask is not None:
        raise SystemExit("Use either --cell-roi or --cell-mask for the whole-cell ROI, not both.")
    if args.model == "diffusion" and args.bleach_mask is not None:
        raise SystemExit("--bleach-mask currently supports reaction workflows; diffusion requires --roi geometry.")

    bleach_roi = CircularROI(*args.roi) if args.roi is not None else None
    if bleach_roi is not None:
        mask_factory = _build_mask_factory(bleach_roi)
        roi_mode = 1
        roi_definition = args.roi
    else:
        mask_factory = _load_mask_argument(args.bleach_mask, "bleach")
        roi_mode = 2
        roi_definition = ()

    cell_factory = None
    if args.cell_mask is not None:
        cell_factory = _load_mask_argument(args.cell_mask, "cell")
    elif args.cell_roi is not None:
        cell_factory = _build_mask_factory(CircularROI(*args.cell_roi))

    adjacent_factory = None
    if args.use_adjacent_roi:
        if bleach_roi is None:
            raise SystemExit("--use-adjacent-roi currently requires circular --roi geometry.")
        adjacent_roi = adjacent_circle(bleach_roi, offset_factor=args.adjacent_offset)
        adjacent_factory = _build_mask_factory(adjacent_roi)

    post_bleach_index = max(args.post_bleach_frame - 1, 0)

    model_index = {"diffusion": 1, "reaction1": 2, "reaction2": 3}[args.model]
    inputs = BasicInputs(
        file_paths=[path.resolve() for path in args.files],
        model_index=model_index,
        roi_mode=roi_mode,
        normalize_by_cell=args.normalize_by_cell,
        background_intensity=args.background,
        post_bleach_frame=post_bleach_index,
        roi_definition=roi_definition,
        pre_bleach_frame_count=args.pre_bleach_count,
        use_adjacent_roi=args.use_adjacent_roi,
    )

    if args.model in {"reaction1", "reaction2"}:
        datasets = load_reaction_datasets(
            inputs,
            bleach_roi=mask_factory,
            cell_roi=cell_factory,
        )
        if not datasets:
            raise SystemExit("No datasets loaded.")

        if args.model == "reaction1":
            config = default_reaction1_config(
                frap_length=len(datasets[0].norm_frap),
                post_bleach_frame=inputs.post_bleach_frame,
            )
            model_order = 1
        else:
            config = default_reaction2_config(
                frap_length=len(datasets[0].norm_frap),
                post_bleach_frame=inputs.post_bleach_frame,
            )
            model_order = 2

        if args.fit_mode == "individual":
            config.fit_averaged_data = False
        elif args.fit_mode == "average_curve":
            config.fit_averaged_data = True

        result = fit_reaction_model(datasets, inputs, config, model_order=model_order)
        _write_cli_exports(args, result, datasets, mask_factory, cell_factory, args.model)

        print(f"Reaction {model_order} model fit results:")
        for name, value in result.parameters.items():
            print(f"  {name}: {value:.4g}")
        print(f"  Photodecay rate: {result.decay_rate:.4g}")
        print(f"  Sum of squared residuals: {result.sum_squared_residuals:.4g}")
        return

    datasets = load_diffusion_datasets(
        inputs,
        bleach_roi=mask_factory,
        adjacent_roi=adjacent_factory,
    )

    if not datasets:
        raise SystemExit("No datasets loaded.")

    config = default_diffusion_config(
        frap_length=len(datasets[0].norm_frap),
        profile_length=len(datasets[0].radius),
        post_bleach_frame=inputs.post_bleach_frame,
    )
    config.fit_mode = args.fit_mode

    result = fit_diffusion_model(datasets, inputs, config)
    _write_cli_exports(args, result, datasets, mask_factory, cell_factory, args.model)

    print("Diffusion model fit results:")
    print(f"  Bleach depth (k): {result.k:.4g}")
    print(f"  Effective radius (re): {result.r_effective:.4g}")
    if result.half_time is not None:
        print(f"  Half time (tau_1/2): {result.half_time:.4g}")
    print(f"  Diffusion coefficient (D): {result.diffusion_coefficient:.4g}")
    print(f"  Mobile fraction (MF): {result.mobile_fraction:.4g}")
    if result.corrected_mobile_fraction is not None:
        print(f"  Corrected mobile fraction: {result.corrected_mobile_fraction:.4g}")
    print(f"  Photodecay rate: {result.decay_rate:.4g}")
    print(f"  Sum of squared residuals: {result.sum_squared_residuals:.4g}")


if __name__ == "__main__":
    main()
