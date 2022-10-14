# Data

## Ensemble Data

IFS data can be downloaded from the ECMWF MARS archive. If you have the required credentials, log in to ecgate and use the `job.sh` script provided in the `IFS-ENS` folder to obtain the files

    $ sbatch job.sh 2016-12-09 12
    $ sbatch job.sh 2016-12-10 00
    $ sbatch job.sh 2016-12-10 12
    $ sbatch job.sh 2016-12-11 00

For the forecast evaluation plot (Figure 5), additional forecasts for 18 Dec 2016 are required.
`IFS-ENS/EVAL` contains the scripts necessary to retrieve this data from the ECMWF archive via ecgate.
Generate a MARS request for each forecast initialization time:

    $ python eval-requests.py

Submit these jobs to the cluster

    $ source eval-submit

When all data has been retrieved, convert the grib files to netcdf.

    $ sbatch eval-to-netcdf

Place the netcdf files into the `IFS-ENS/EVAL` folder.

Make sure to specify appropriate paths in the scripts first.

Data may alternatively be downloaded via the [ECMWF API client](https://github.com/ecmwf/ecmwf-api-client) after translation of the MARS request to a web API request (not provided). If you don't have access rights to operational ensemble forecast data, [TIGGE](https://confluence.ecmwf.int/display/TIGGE/TIGGE+archive) offers similar data at lower vertical resolution.


## Reanalysis Data

ERA5 data can be downloaded from the Copernicus Climate Data Store with the provided script

    $ python3 download-ERA5.py`

(call from within this directory).

