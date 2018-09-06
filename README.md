# "Swimming downstream" evaluation code

This repository contains all the code used in the evaluation of
methods for DTU, DGE, and DTE in the F1000Research article,
*Swimming downstream: statistical analysis of differential transcript*
*usage following Salmon quantification*. The directories contain the
following analyses:

* `countsimqc` - generate the countsimQC report (and the report)
* `dte` - generate the DGE and DTE evaluations
* `dtu` - generate the DTU evaluation
* `quant` - bash scripts to quantify the abundances 
* `rats` - run RATs
* `simulate` - generate the simulated reads
* `suppa` - run SUPPA2

Note on the code in `simulate`: simulating reads proceeds in three
steps: 1) `simulate_expression.R` is used to simulate expression
vectors for two groups with DGE, DTE, and DTU genes, 2)
`simulate_reads.R` is used to actually simulate the paired-end reads
with the *polyester* package, 3) `shuffle_wrapper.sh` is used to
shuffle and then compress the paired-end reads output by
*polyester*. The simulated reads need to be shuffled to be input to
*Salmon* (because *polyester* outputs reads from each transcript, one
at a time).

*DRIMSeq* and *DEXSeq* were run using the code directly from the
`rnaseqDTU` workflow, see <https://github.com/mikelove/rnaseqDTU>.
