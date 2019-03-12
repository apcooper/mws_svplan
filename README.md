# MWS SV planning notebooks 

These notebooks are experiments with the tile layout, target selection, priority scheme etc. for MWS SV. For a given plan, each experiment assigns fibers to targets.

Changing target selections and priorities requires some hacking of `desitarget` (see Background Information below).
 
## Summary

Experiments with different target selections are split into pairs of notebooks: `sv_design_$RUN.ipynb` and `sv_analysis_$RUN.ipynb` where `$RUN` is a unique label for the experiment and corresponds to an output dir under `./runs`.

The 'design' notebook deals with the tile setup and the generation of data. The 'analysis' notebook is for the analysis of the data. The split into two notebooks is so that all the analysis steps can be re-run without regenerating the data.

The directories in this repository are:

* `./cache`: local views of large data files (e.g. sweeps).
* `./py`: python module with various helper functions.
* `./runs`: data output by notebooks goes here.
* `./scripts`: code fragments for use in notebooks.
* `./targetmasks`: custom versions of desitarget.yaml.

## Instructions

Each notebook starts with a cell that imports a preamble common to all experiments: 
```python
# %load scripts/nbpreamble.py
...
```

Uncomment this line and run the cell. This only needs to be done if the contents of `scripts/nbpreamble` have changed, otherwise the cell can just be run as-is. Sometimes this cell needs to be run twice to get it to work if the preamble is updated.

## Tutorials 

## Background Information / caveats

The `desitarget` code defines the MWS target selction as a set of python fuctions. These functions need to be adapted for SV.

The approach here is to write drop-in replacements for the relevant functions (done in `py/apcsv/cuts.py`) and then patch them by hand into whatever version of `desitarget` is imported by default for your chosen notebook kernel. 

This approach is a bit abstruse. The alternative is to have a separate branch of the full `desitarget` code and import that instead of the default version. For casual use, the patching approach has less user setup overhead (in getting jupyter to import the right version and pick up changes) and works well enough, so long as the patched functions also track any unrelated changes to the same functions in desitarget (so don't use this with old versions of desitarget).

Since the target bitmask is read from `desitarget/py/data` on import, the notebooks also need to patch that with a custom version. See `patch_desitarget_bits()` in `py/apcsv/patch.py`.


