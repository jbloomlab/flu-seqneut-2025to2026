# README
This repository designs HA insert sequences that can be submitted to Twist Biosciences for gene fragment synthesis and cloning. These inserts are designed to fit between the BsmBI-v2 and XbaI cut sites in the Bloom lab vector 2851. They should begin after the 19th codon of WSN upstream signal peptide, and continue all the way through the end of the HA coding region, followed by a double stop codon and a 16-nucleotide barcode. 

In our design, our H3 HA constructs have an upstream signal peptide from WSN, an HA ectodomain matching a currently circulating strain, an endodomain matching a recent H3 consensus sequence, and a C-terminal domain from WSN

Our H1 HA constructs are slightly different, with an upstream signal peptide from WSN, an HA ectodomain matching a currently circulating strain, and an endodomain and a C-terminal domain also from WSN

We specifically design our barcodes to avoid barcodes used by prior libraries and barcodes starting with `GG` due to sequencing issues. Prior libraries are identified in the `configuration.yml`. 

To run the pipeline, build the environment with:
```
conda env create -f environment.yml
```

Then run the analysis notebook `construct_order.py` interactively in Marimo. 