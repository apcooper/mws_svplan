# MWS SV planning notebooks 

These notebooks are experiments with the tile layout, target selection, priority scheme etc. for MWS SV. For a given plan, each experiment assigns fibers to targets.

Changing target selections and priorities requires some hacking of `desitarget` (see Background Information).

## Instructions

* `/cache`: local views of large data files (e.g. sweeps).
* `/py`: python module with various helper functions.
* `/runs`: data output by notebooks goes here.
* `/scripts`: code fragments for use in notebooks.
* `/targetmasks`: custom versions of desitarget.yaml.

Experiments are split into pairs of notebooks: `sv_design_$RUN.ipynb` and `sv_analysis_$RUN.ipynb` where `$RUN` is a unique label for the experiment and corresponds to an output dir under `/runs`.

The 'design' notebook deals with the tile setup and the generation of data. The 'analysis' notebook is for the analysis of the data.

## Background Information

The `desitarget` code defines the MWS selction function using several python fuctions. The approach in these notebook is to write drop-in replacements for those functions (in `py/apcsv/cuts.py`) and patch them in to whatever version of `desitarget` is imported by default for your notebook kernel. This is a bit abstruse. The (better) alternative is to have a separate branch of the full `desitarget` code and import that instead of the default version.

Since the target bitmask is read from `desitarget/py/data` on import, the notebooks also need to patch that with a custom version. See `patch_desitarget_bits()` in `py/apcsv/patch.py`.


