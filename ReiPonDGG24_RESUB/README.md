Code, Data, and Results for Numerical Experiments in "$\mathcal{H}_2$ 
optimal model reduction of linear systems with multiple quadratic outputs"
===========================================================================

This archive contains the companion codes, data and computed results for 
the paper:

S. Reiter, I. Pontes Duff, I. V. Gosea, S. Gugercin; "$\mathcal{H}_2$ 
optimal model reduction of linear systems with multiple quadratic outputs",

which implement numerical experiments using $\mathcal{H}_2$ model order 
reduction methods for linear quadratic-output systems.


## Dependencies and Installation

The code was tested under MATLAB 2023b.
The Matrix Equation Sparse Solver (M-M.E.S.S.) Library version 3.0. is 
used for solving the sparse-dense Sylvester equations in 
the `\drivers\mimolqo_tsia.m` function.


## Getting Started

The file `runme_advecdiff.m` can be used to reproduce experiments from the 
corresponding paper on the advection diffusion problem with a quadratic 
cost function as the quantity of interest. The model can be found under
`data\AdvecDiff_n3000.mat`. The file `runme_advecdiff.m` generates the 
reduced-order models computed using `\drivers\mimolqo_tsia.m` and 
`\drivers\mimolqo_bt.m`, and reproduces the experiments found in Section 5 
of the companion paper.

The results computed by these scripts will be saved to the `results`
folder. Existing results will be overwritten.


## Author

Sean Reiter
* affiliation: Virginia Tech (USA)
* email: seanr7@vt.edu
* orcid: [0000-0002-7510-1530](https://orcid.org/0000-0002-7510-1530)


## License

Copyright (C) 2024 Sean Reiter

In general, this software is licensed under the BSD-2 License.
See [COPYING](COPYING) for a copy of the license.

The files in `results` are licensed under the CC-BY 4.0 License.
See [COPYING_DATA](COPYING_DATA) for a copy of the license.


## Citation


### DOI

The DOI for version 1.2 is
[10.5281/zenodo.14532348](https://doi.org/10.5281/zenodo.14532348).


### Cite as

S. Reiter. Code, data, and results for numerical 
experiments in "$\mathcal{H}_2$ optimal model reduction of linear systems 
with multiple quadratic outputs" (version 1.2),
December 2024. doi:10.5281/zenodo.11104814


### BibTeX

    @MISC{supRei24,
      author =       {Reiter, Sean},
      title  =       {Code, Data, and Results for Numerical Experiments in 
                      ``$\mathcal{H}_2$ optimal model reduction of linear 
                      systems with multiple quadratic outputs'' 
                      (version 1.2)},
      month  =       dec,
      year   =       {2024},
      doi    =       {10.5281/zenodo.14532348}
    }