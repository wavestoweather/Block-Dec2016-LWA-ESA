# Block-Dec2016-LWA-ESA

Polster, C., and V. Wirth, 2023: The Onset of a Blocking Event as a “Traffic Jam”: Characterization with Ensemble Sensitivity Analysis. *J. Atmos. Sci.*, https://doi.org/10.1175/JAS-D-21-0312.1, in press.

This repository contains the code to produce all figures of the article.

- Research by [Christopher Polster](https://dynmet.ipa.uni-mainz.de/christopher-polster/) and [Volkmar Wirth](https://dynmet.ipa.uni-mainz.de/volkmar-wirth/).
- Software by Christopher Polster.
- [hn2016_falwa](https://github.com/csyhuang/hn2016_falwa/) package by [Clare S. Y. Huang](https://csyhuang.github.io/). The version included in this repository is 0.5.0 from August 2021. The package has seen significant improvements since (see [releases](https://github.com/csyhuang/hn2016_falwa/releases)) and we recommend using the latest version.

[LICENSE](LICENSE). The contents of `scripts/hn2016_falwa` are re-distributed under the terms of their [LICENSE](scripts/hn2016_falwa/LICENSE.txt).

[![DOI](https://zenodo.org/badge/430733123.svg)](https://zenodo.org/badge/latestdoi/430733123)


## How To Run

Clone this repository:

    $ git clone https://github.com/wavestoweather/Block-Dec2016-LWA-ESA.git

### Data Preparation

Use the scripts and instructions in the [`data`](data) directory to obtain the required input data. Before continuing, the `data` folder should look like this:

    data
    + IFS-ENS
    | + ENS-2016-12-09T12Z
    | | + ENS-2016-12-09T12Z-t.nc
    | | + ENS-2016-12-09T12Z-u.nc
    | | + ENS-2016-12-09T12Z-v.nc
    | + ENS-2016-12-10T00Z
    | | + ENS-2016-12-10T00Z-t.nc
    | | + ENS-2016-12-10T00Z-u.nc
    | | + ENS-2016-12-10T00Z-v.nc
    | + ENS-2016-12-10T12Z
    | | + ENS-2016-12-10T12Z-t.nc
    | | + ENS-2016-12-10T12Z-u.nc
    | | + ENS-2016-12-10T12Z-v.nc
    | + ENS-2016-12-11T00Z
    | | + ENS-2016-12-11T00Z-t.nc
    | | + ENS-2016-12-11T00Z-u.nc
    | | + ENS-2016-12-11T00Z-v.nc
    | + EVAL
    | | + ENS-DEC18-EVAL-2016-12-08T00Z.nc
    | | + ENS-DEC18-EVAL-2016-12-08T12Z.nc
    | | + ...
    | | + ENS-DEC18-EVAL-2016-12-18T00Z.nc
    | | + eval-requests.py
    | | + eval-submit
    | | + eval-to-netcdf
    | + job.sh
    + download-ERA5.py
    + ERA-2016-12-10-to-2016-12-24-uvt.nc
    + README.md

Note that while the analysis and plots are based on 2° data, 1° data is obtained by the provided download scripts.
Input data is coarsened during processing (`--half-resolution` parameter of `scripts.calculate-lwa`).

### Software Requirements

The following software needs to be available to run the data processing and plotting scripts:

- make
- Fortran compiler (for building the `hn2016_falwa` Python extensions)
- C compiler with OpenMP (for building Python extensions)
- Python 3 with packages listed in [`requirements.txt`](requirements.txt) (the specified minumum versions were used during development, older versions may or may not work)
- lualatex with fontspec, amsmath, tikz, graphicx (Fig. 1 and to combine panels for Fig. 8)

### Data Analysis and Plots

Run

    $ make

from the top-level of the repository. This will compile the Python extensions, compute the intermediate PV, LWA and flux fields and create all plots. Use `--dry-run` to check the sequence of commands to be executed by `make` first. Use `-jN` to run `N` jobs in parallel if you have sufficient CPU and memory resources to do so.

### Output

Figures will appear in the `figures` directory:

- Fig. 1: `figures/forecast-combination.pdf`
- Fig. 2: `figures/event24-reanalysis+nh18fig5.pdf`
- Fig. 3: `figures/event24-plume.pdf`
- Fig. 4: `figures/event24-budget.pdf`
- Fig. 5: `figures/event24-evaluation.pdf`
- Fig. 6: `figures/idealized.pdf`
- Fig. 7: `figures/event24-maps+hovmoeller.pdf`
- Fig. 8: `figures/event24-cluster+scatter.pdf`
- Fig. 9: `figures/event24-maps-separate.pdf`
- Fig. 10: `figures/event24-cluster-f2.pdf`
- Fig. 11: `figures/event24-nh18fig4.pdf`, does not contain the fits of [Nakamura and Huang (2018)](https://doi.org/10.1126/science.aat0721) due to licensing restrictions

Because the statistical significance test is based on a randomized procedure, results may vary slightly each time the plots are created.


## Acknowledgements

This research project has been carried out within the Transregional Collaborative Research Center [SFB/TRR 165](https://www.dfg.de/en/funded_projects/current_projects_programmes/list/projectdetails/index.jsp?id=257899354) "Waves to Weather" funded by the German Science Foundation (DFG). https://www.wavestoweather.de/

