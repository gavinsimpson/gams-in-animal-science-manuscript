# 📝 Using generalized additive models in the livestock animal sciences

## License

The R code in the repo is released under an MIT license. The manuscript (`manuscript.qmd` and `manuscript.pdf`) and preprint (`preprint-version.qmd` and `preprint-version.pdf`) and the figures in the `manuscript_files/` and `preprint-version_files/` folders are made available under a the [Creative Commons By-Attribution v4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode.en) license (CC-BY).

## Description

This repo contains the Quarto sources for my manuscript "Method: Using generalized additive models in the livestock animal sciences", which has been submitted to Animal &mdash; open space. The aim of this work is to promote the understanding and adoption of GAMs in animal science, where nonlinear covariate–response relationships are common.

The revised version responds to reviewer feedback requesting a clearer demonstration of why researchers in animal science should consider generalized additive models (GAMs) and whether their use extends beyond prediction.

To address this, I have extended the examples to illustrate how posterior simulation from fitted GAMs can be used for biologically meaningful inference:

* 🐄 Using the modelled lactation curve, I now show how to estimate the timing and magnitude of peak milk fat production, and the persistency of milk fat content in late lactation, and compare the estimates with those of more traditional methods.

* 🐖 Deriving growth rate estimates for pigs from hierarchical GAMs fitted to depth-camera weight data. For this we need the derivative of the fitted growth curves on the response (pig weight) scale and posterior sampling enables easy estimation of the uncertainties in the growth rates.

* 🪶 This is in addition to the formal pairwise comparisons of quail weight and growth rates in an experiment on maternal thyroid hormones, which were already in the manuscript, using the excellent [marginaleffects](https://marginaleffects.com/) 📦 for R

These examples demonstrate how GAMs can provide flexible, data-driven inference that yields interpretable quantities of direct biological relevance.

The `analysis` folder contains three R scripts which work through the examples in the paper.

The `manuscript.qmd` file is the version of the manuscript submitted to Animal &mdash; open space, while `preprint-version.qmd` is a version that was submitted to the ArXiv as a preprint (the preprint has figures and tables within the text, not at the end, and I fixed a typo or two while I was making those changes, but in all other respects the two versions are the same).
