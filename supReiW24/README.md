Code, Data and Results for Numerical Experiments in "Interpolatory 
model order reduction of large-scale dynamical systems with root mean 
squared error measures"
===========================================================================

This archive contains the companion codes, used data and computed results
for the paper:

S. Reiter, S. W. R. Werner; "Interpolatory model order reduction of 
large-scale dynamical systems with root mean squared error measures",

which implement numerical experiments with different interpolatory model 
order reduction methods for linear quadratic output systems that 
model root mean squared error measures.


## Dependencies and Installation

The code was tested under MATLAB 2023b.


## Getting Started

The `runme*.m` files can be used to reproduce experiments from the
corresponding paper.
The scripts correspond to the following experiments:
* `runme`: experiments using the vibro-acoustic plate with tuned vibration 
  absorbers (plateTVA) data set
* `runme_benchmarks`: code to generate the interpolatory reduced-order 
  models computed using `\drivers\interpolatory_solves.m` and used in 
  `runme.m`
* `runme_fosim`: code to run the full-order simulation of the plateTVA data
  set
* `runme_lqoirka`: code to generate the interpolatory reduced-order models 
  computed using `\drivers\sisolqo_irka.m` and used in `runme.m`

The results computed by these scripts will be saved to the `results`
folder. Existing results will be overwritten.

The plateTVA data is licensed to Quirin Aumann under the CC-BY 4.0 License 
and available at https://zenodo.org/records/7671686


## Author

Sean Reiter
* affiliation: Virginia Tech (USA)
* email: seanr7@vt.edu
* orcid: [0000-0002-7510-1530](https://orcid.org/0000-0002-7510-1530)

Steffen W. R. Werner
* affiliation: Virginia Tech (USA)
* email: steffen.werner@vt.edu
* orcid: [0000-0003-1667-4862](https://orcid.org/0000-0003-1667-4862)


## License

Copyright (C) 2024 Sean Reiter, Steffen W. R. Werner

In general, this software is licensed under the BSD-2 License.
See [COPYING](COPYING) for a copy of the license.

The files in `results` and `data` with its content are licensed under 
the CC-BY 4.0 License. See [COPYING_DATA](COPYING_DATA) for a copy of the 
license.


## Citation


### DOI

The DOI for version 1.0 is
[10.5281/zenodo.10729524](https://doi.org/10.5281/zenodo.10729524).


### Cite as

S. Reiter and S. W. R. Werner. Code, Data and Results for Numerical 
Experiments in "Interpolatory model order reduction of large-scale 
dynamical systems with root mean squared error measures" (version 1.0),
March 2024. doi:10.5281/zenodo.1072952


### BibTeX

    @MISC{supReiW24,
      author =       {Reiter, S. and Werner, S.~W.~R.},
      title  =       {Code, Data and Results for Numerical 
                      Experiments in ``{I}nterpolatory model order 
                      reduction of large-scale dynamical systems with root 
                      mean squared error measures'' (version 1.0)},
      month  =       mar,
      year   =       {2024},
      doi    =       {10.5281/zenodo.zenodo.1072952}
    }
