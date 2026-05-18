from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

from .data.loading import load_diffusion_datasets
from .models.diffusion import default_diffusion_config, fit_diffusion_model
from .roi import CircularROI, adjacent_circle
from .types import BasicInputs


def _build_mask_factory(roi: CircularROI):
    def factory(shape):
        return roi.to_mask(shape)

    return factory


def parse_arguments(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="FRAP-Toolbox Python port")
    parser.add_argument("files", nargs="+", type=Path, help="Paths to FRAP image files")
    parser.add_argument("--model", choices=["diffusion"], default="diffusion")
    parser.add_argument(
        "--roi",
        nargs=3,
        type=float,
        metavar=("X", "Y", "RADIUS"),
        required=True,
        help="Circular ROI definition in pixels",
    )
    parser.add_argument("--post-bleach-frame", type=int, default=21)
    parser.add_argument("--pre-bleach-count", type=int, default=10)
    parser.add_argument("--background", type=float, default=0.0)
    parser.add_argument("--normalize-by-cell", action="store_true")
    parser.add_argument("--use-adjacent-roi", action="store_true")
    parser.add_argument(
        "--adjacent-offset",
        type=float,
        default=2.5,
        help="Offset multiplier for adjacent circular ROI (used with --use-adjacent-roi)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_arguments(argv)

    bleach_roi = CircularROI(*args.roi)
    mask_factory = _build_mask_factory(bleach_roi)

    adjacent_factory = None
    if args.use_adjacent_roi:
        adjacent_roi = adjacent_circle(bleach_roi, offset_factor=args.adjacent_offset)
        adjacent_factory = _build_mask_factory(adjacent_roi)

    post_bleach_index = max(args.post_bleach_frame - 1, 0)

    inputs = BasicInputs(
        file_paths=[path.resolve() for path in args.files],
        model_index=1,
        roi_mode=1,
        normalize_by_cell=args.normalize_by_cell,
        background_intensity=args.background,
        post_bleach_frame=post_bleach_index,
        roi_definition=args.roi,
        pre_bleach_frame_count=args.pre_bleach_count,
        use_adjacent_roi=args.use_adjacent_roi,
    )

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

    result = fit_diffusion_model(datasets, inputs, config)

    print("Diffusion model fit results:")
    print(f"  Bleach depth (k): {result.k:.4g}")
    print(f"  Effective radius (re): {result.r_effective:.4g}")
    print(f"  Diffusion coefficient (D): {result.diffusion_coefficient:.4g}")
    print(f"  Mobile fraction (MF): {result.mobile_fraction:.4g}")
    if result.corrected_mobile_fraction is not None:
        print(f"  Corrected mobile fraction: {result.corrected_mobile_fraction:.4g}")
    print(f"  Photodecay rate: {result.decay_rate:.4g}")
    print(f"  Sum of squared residuals: {result.sum_squared_residuals:.4g}")


if __name__ == "__main__":
    main()
