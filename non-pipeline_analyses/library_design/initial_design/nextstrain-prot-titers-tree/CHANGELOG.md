# CHANGELOG

### version 1.1.0
+ Build the trees faster with `iqtree`
  - The tree building uses multiple threads if specified when running snakemake.
  - Use the `-fast` option to `iqtree`

+ Handle if `treetime` cannot infer a date

## version 1.0.0
Initial version.
