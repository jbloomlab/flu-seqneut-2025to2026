# Build interactive `nextstrain` trees on protein sequences designed to display neutralization titer values

[![Release](https://img.shields.io/github/v/release/jbloomlab/nextstrain-prot-titers-tree?logo=github)](https://github.com/jbloomlab/nextstrain-prot-titers-tree/releases)
[![Build Status](https://github.com/jbloomlab/nextstrain-prot-titers-tree/actions/workflows/test.yaml/badge.svg)](https://github.com/jbloomlab/nextstrain-prot-titers-tree/actions/workflows/test.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

---

This repository contains a `snakemake` pipeline developed by the [Bloom lab](https://jbloomlab.org) that builds interactive `nextstrain` trees of protein sequences that can be colored and analyzed in terms of additional data such as neutralization titers.
The pipeline was designed for the use case of displaying high-throughput neutralization titer data for many strains similar to that described in [Kikawa et al (2025)](https://doi.org/10.1101/2025.09.06.674661).

This pipeline is specifically tailored for the case where you want to build **protein** sequence trees and have the divergence indicate the number of amino-acid mutations separating different proteins.
More standard `nextstrain augur` pipelines may be more appropriate if you are using nucleotide sequences.

## Configuring the pipeline, running it, and viewing the results
To run the pipeline, you need to build a configuration pipeline that has the configuration for the tree (input data, display options, etc).
See [example_config.yaml](example_config.yaml), which has an example configuration using the H3N2 data from [Kikawa et al (2025)](https://doi.org/10.1101/2025.09.06.674661) as stored in [./example_data/](example_data).
You should build your own configuration file for your data mirroring that example.
Then run the pipeline with:

    snakemake -j <nthreads> --configfile <path_to_your_configuration_file> --software-deployment-method conda

Note that running this requires `snakemake` to be installed, which you can do by building and activating the `conda` environment in [environment.yml](environment.yml).

The tree-building step using IQ-TREE will use multiple threads (up to a maximum of 8 threads, or the number of cores specified with the `-j` argument to `snakemake`, whichever is smaller) to speed up the analysis.

The result of this is an auspice JSON file with the tree suitable for viewing either by uploading to [https://auspice.us/](https://auspice.us/) or via a [Nextstrain Community Build](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html).
The auspice JSON tree for the example is in [auspice/nextstrain-prot-titers-tree.json](auspice/nextstrain-prot-titers-tree.json) and can be viewed as a [Nextstrain Community Build](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html) at [https://nextstrain.org/community/jbloomlab/nextstrain-prot-titers-tree@main](https://nextstrain.org/community/jbloomlab/nextstrain-prot-titers-tree@main).

If the *metadata* in the configuration file has titers, they are displayed on the tree as in the above example.
You can also show all amino-acid identities on the tree, color by amino-acid identity at a site, and show branch lengths either based on amino-acid mutations per site or time.

If you also specify *titers* with per-serum titers (as indicated in [example_config.yaml](example_config.yaml)) then the pipeline will also produce a sidecar JSON with these titers (eg, for this example [auspice/nextstrain-prot-titers-tree_measurements.json](auspice/nextstrain-prot-titers-tree_measurements.json)) that can be used to analyze per-serum titers in the *Measurements* panel when viewing the tree.

## Using in a larger `snakemake` pipeline
The typical way to use this pipeline is as a submodule of a larger `snakemake` pipeline.
See [https://github.com/jbloomlab/flu-seqneut-2025](https://github.com/jbloomlab/flu-seqneut-2025) for an example of how that can be done.

Briefly, first add this repo as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to your larger repository pipeline by cloning it into that repository and then additing it as a git submodule with:

      git submodule add https://github.com/jbloomlab/nextstrain-prot-titers-tree

This creates a file called `gitmodules` and adds the `nextstrain-prot-titers-tree` subdirectory, both of which can then be committed to your parent repo.

You can then use it as a module in your larger pipeline, as for instance like this:

```
for subtype in config["subtypes"]:
    module_name = f"nextstrain-prot-titers-tree_{subtype}"
    module:
        name: module_name
        snakefile: "nextstrain-prot-titers-tree/Snakefile"
        config: config["nextstrain-prot-titers-tree_config"][subtype]
    use rule * from module_name as module_name*
```


## Testing via GitHub Actions
When updating the pipeline, you should:

 - lint code with [ruff](https://github.com/astral-sh/ruff) (`ruff check .`)
 - format code with [black](https://github.com/psf/black) (`black .`)
 - lint [Snakefile](Snakefile) with [snakemake --lint](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html) (`snakemake --lint --configfile example_config.yaml`)
 - format [Snakefile](Snakefile) with [snakefmt](https://github.com/snakemake/snakefmt) (`snakefmt .`).

These checks are run automatically when you via the GitHub Action specified in [.github/workflows/test.yaml](.github/workflows/test.yaml).
