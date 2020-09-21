# Structural Eigenmodes of Complex Graph Laplacian
---
This repository includes the analysis as described in [Emergence of Canonical Functional Networks from Complex-valued Laplacian of Structural Connectome](https://www.biorxiv.org/content/10.1101/2020.09.16.300384v1), where we investigate whether we can recreate the canonical functional networks of the human brain from structural connectivity with a simple low dimensional framework. 

## Set up:
This repository is dependent on [`spectrome`](https://github.com/Raj-Lab-UCSF/spectrome), please see instructions in the repository for setting up your `conda` environment. The `spectrome` repository will need to be cloned, and the path to the `spectrome` folder will need to be appended to `$PYTHONPATH`. Activating the `spectrome` environment will enable all the analysis in this repository.

If you wish to visualize brain renderings, you will need to add [`ipywidgets`](https://ipywidgets.readthedocs.io/en/latest/) and [`pysurfer`](https://pysurfer.github.io/#more-information) dependencies to your environment. The current public `spectrome` environment is made for Binder, meaning we removed any pop out visualization tools that caused issues with Binder. 

## Citation:


## Files:
 - `notebooks`: Jupyter notebooks that produced the paper results.
 - `scripts`: python script for optimization. 
 - `data`: some intermediate data needed to generate the figures.
    - `mean80_fiber*.csv`: HCP averaged fibercount and fiberlength files. 
    - `DK_dictionary_normalized.csv`: The canonical functional networks in DK-atlas parcellations