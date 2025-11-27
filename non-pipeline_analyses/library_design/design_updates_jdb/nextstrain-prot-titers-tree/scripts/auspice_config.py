import json
import sys

import matplotlib

import numpy

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# get discrete colors spaced across color map
n_scale_points = 8  # use this many points to define the color scale
hex_colors = {
    cmap: [
        matplotlib.colors.to_hex(matplotlib.cm.get_cmap(cmap, n_scale_points)(i))
        for i in range(n_scale_points)
    ]
    for cmap in ["viridis", "viridis_r"]
}


metadata = pd.read_csv(snakemake.input.metadata, sep="\t")
color_by_metadata = snakemake.params.color_by_metadata
assert set(color_by_metadata).issubset(metadata.columns)

colorings = []
for col, col_d in color_by_metadata.items():
    if "exclude_auto_scale" in col_d:
        assert isinstance(col_d["exclude_auto_scale"], list), col_d[
            "exclude_auto_scale"
        ]
        for_lims = metadata[~metadata["strain"].isin(col_d["exclude_auto_scale"])]
    else:
        for_lims = metadata
    minval = col_d["fixed_min"] if ("fixed_min" in col_d) else for_lims[col].min()
    maxval = col_d["fixed_max"] if ("fixed_max" in col_d) else for_lims[col].max()
    if minval > metadata[col].min():
        minprefix = "<="
    else:
        minprefix = ""
    if maxval < metadata[col].max():
        maxprefix = ">="
    else:
        maxprefix = ""
    assert maxval >= minval, f"{maxval=}, {minval=}"
    if "scale_type" in col_d:
        scale_type = col_d["scale_type"]
        if scale_type.endswith(("_linear", "_log")):
            scale = scale_type.split("_")[-1]
            cmap = "_".join(scale_type.split("_")[:-1])
        else:
            raise ValueError(f"invalid {scale_type=}")
        if cmap in hex_colors:
            hexcols = hex_colors[cmap]
        else:
            raise ValueError(f"invalid {scale_type=}")
        if scale == "linear":
            scalevals = numpy.linspace(minval, maxval, num=len(hexcols))
        elif scale == "log":
            scalevals = numpy.logspace(
                numpy.log(minval) / numpy.log(2),
                numpy.log(maxval) / numpy.log(2),
                num=len(hexcols),
                base=2,
            )
        else:
            raise ValueError(f"{scale=}")
        assert scale_type in [
            "viridis_linear",
            "viridis_log",
            "viridis_r_linear",
            "viridis_r_log",
        ]
        legendlabels = [f"{v:.3g}" for v in scalevals]
        legendlabels[0] = minprefix + legendlabels[0]
        legendlabels[-1] = maxprefix + legendlabels[-1]
        color_scale = {
            "type": "continuous",
            "scale": [[float(v), c] for v, c in zip(scalevals, hexcols)],
            "legend": [
                {"value": v, "display": d} for v, d in zip(scalevals, legendlabels)
            ],
        }
    else:
        color_scale = {}
    colorings.append({"key": col} | color_scale)

json_d = {
    "display_defaults": snakemake.params.display_defaults,
    "colorings": colorings,
    "panels": (
        ["tree", "entropy"] + (["measurements"] if snakemake.params.has_titers else [])
    ),
}

with open(snakemake.output.auspice_config, "w") as f:
    json.dump(json_d, f, indent=2)
