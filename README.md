# MTX_synthetic User Manual #

Thid module contains Python code for generating synthetic metatranscriptomic (MTX) and metagenomic (MGX) count data based on template microbial community taxonomic abundances. If you use the MTX_synthetic code, please cite our manuscript:

Zhang et al., "Statistical approaches for differential expression analysis in metatranscriptomics" (In submission).

And feel free to link it to your Methods: https://github.com/biobakery/MTX_synthetic.

Contact Eric Franzosa (eric.franzosa@gmail.com) for help.

# Requirements

- Python 2.7+ or Python 3+
- Python `numpy`
- Python `doit`
  - https://pydoit.org/
  - only needed to reproduce existing datasets
- Python `zopy` 
  - https://github.com/franzosa/zopy
  - Eric Franzosa's Python codebase
  - Contact eric.franzosa@gmail.com for help

# Inputs

- Unzip the included HMP1-II MetaPhlAn 2 taxonomic abundances

# Usage

- Execute `python synth_mgx_mtx.py <hmp> <site> <basename>`
  - `<hmp>` is the path to the abundance file
  - `<site>` is an HMP environment type, e.g. "Stool"
  - `<basename>` is a prefix for the files to be generated
- Consult the script's help menu for tuning the parameters of the simulation

```
usage: synth_mgx_mtx.py [-h] [--n-bugs N_BUGS] [--n-samples N_SAMPLES]
                        [--n-groups N_GROUPS] [--n-pangenes N_PANGENES]
                        [--coreness CORENESS]
                        [--dna-read-depth DNA_READ_DEPTH DNA_READ_DEPTH]
                        [--rna-read-depth RNA_READ_DEPTH RNA_READ_DEPTH]
                        [--spike-dep]
                        [--spike-dep-strength SPIKE_DEP_STRENGTH]
                        [--spike-groups] [--spike-bug SPIKE_BUG]
                        [--spike-bug-strength SPIKE_BUG_STRENGTH]
                        [--spike-enc SPIKE_ENC]
                        [--spike-enc-strength SPIKE_ENC_STRENGTH]
                        [--spike-exp SPIKE_EXP]
                        [--spike-exp-strength SPIKE_EXP_STRENGTH]
                        hmp site basename
                        ```

# Outputs

A given run produces the following files:

- `<basename>.bug_abunds.tsv`
  - contains species x sample MGX count data
- `<basename>.mgx_abunds.tsv
  - contains gene x sample MGX count data 
- `<basename>.mgx_abunds_groups.tsv
  - contains gene x sample MGX count data, summed over gene families 
- `<basename>.mtx_abunds.tsv
  - contains gene x sample MTX count data
- `<basename>.mtx_abunds_groups.tsv
  - contains gene x sample MTX count data, summed over gene families
- `<basename>.*_spiked.tsv
  - zero or more lists of spiked features

# Doit Workflow

The included `dodo.py` workflow file contains sample script calls for building 11 synthetic meta-omic datasets based on HMP Stool with optional species, gene, transcript, and sequencing depth spikes against a synthetic case:control phenotype (as analyzed in the publication referenced above). Execute this workflow with the `doit` command to build these datasets and their component files.