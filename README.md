# MWS SV planning notebooks 

These notebooks are experiments with the tile layout, target selection, priority scheme etc. for MWS SV. For a given plan, each experiment assigns fibers to targets.

Changing target selections and priorities requires some hacking of `desitarget` (see Background Information below). The aim of the setup here is to make that as painless as possible, so that we can quickly make trials with different selection functions. T
 
## Summary

Experiments with different target selections are split into pairs of notebooks in the toplevel directory of the repository: `$RUN_design.ipynb` and `$RUN_analysis.ipynb` where `$RUN` is a unique label for the experiment and corresponds to an output dir under `./runs`.

The 'design' notebook deals with the tile setup and the generation of data. The 'analysis' notebook is for the analysis of the output of fiberassignment. See examples below. The split into two notebooks is so that all the analysis steps can be re-run without regenerating the data.

The directories in this repository are:

* `./cache`: local subsets of large data files (e.g. sweeps).
* `./py/apc_sv`: python module with various helper functions.
* `./runs`: data output by notebooks goes here.
* `./scripts`: code fragments for use in notebooks.
* `./targetmasks`: custom versions of desitarget.yaml.

## Instructions

All notebooks should use one of the desiconda kernels (18.2, master, etc.).

### Per-notebook setup

Each notebook starts with a cell that imports a preamble common to all experiments: 
```python
# %load scripts/nbpreamble.py
...
```

Uncomment this line and run the cell. This only needs to be done if the contents of `scripts/nbpreamble` have changed, otherwise the cell can just be run as-is. Sometimes this cell needs to be run twice to get it to work if the preamble is updated.

## Examples

### `example_mainsurvey_design.ipynb`

The `design` notebook does the following:

* Makes a cache of sweep headers for a particular set of sweeps. Set `SWEEPS_RESET = True` on first run of (any) notebook and `false` thereafter, only need to change again if using different DR for sweeps. Different notebokes share the cahce of sweep hearders.
* Defines a set of DESI pointings to observe (not bound to the standard tile set) -- functions in the `apc_sv` python package do a quick setup of single-FOV pointings and simple dithers;
* Finds the sweep bricks that overlap those tiles;
* Creates a target catalog from those sweep bricks using `select_targets`:
  - Patches the `targetmask` target definitions and the selection functions called by `select_targets`
* Creates an MTL file (fiberassign input) from the target catalog;
* Selects standard stars from the MTL;
* Runs `fiberassign`.

In this example, I use the main survey selection function and include BGS targets, so the output is representative of one tile from pass 1 of the main survey. 

### `example_mainsurvey_analysis.ipynb`

The intention is to do some useful analysis of the output data from fiberassign. Currently this just makes a few simple plots to demonstrate setting up the notebook and parsing the output.

## Background Information / caveats

The `desitarget` code defines the MWS target selction as a set of python fuctions. These functions need to be adapted for SV.

The approach here is to write drop-in replacements for the relevant functions (done in `py/apcsv/cuts.py`) and then patch them by hand into whatever version of `desitarget` is imported by default for the desiconda notebook kernel you're using. 

This approach is a bit abstruse. The alternative is to have a separate branch of the full `desitarget` code and import that instead of the default version. For casual use, the patching approach has less user setup overhead (in getting jupyter to import the right version and pick up changes) and works well enough, so long as the patched functions also track any unrelated changes to the same functions in desitarget (so don't use this with old versions of desitarget).

Since the target bitmask is read from `desitarget/py/data` on import, the notebooks also need to patch that with a custom version. See `patch_desitarget_bits()` in `py/apcsv/patch.py`.


