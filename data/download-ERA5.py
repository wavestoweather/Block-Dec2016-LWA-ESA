import cdsapi
import datetime as dt


start = dt.datetime(2016, 12, 10)
end   = dt.datetime(2016, 12, 24)


if __name__ == "__main__":

    # List of dates to fetch
    dates = [start]
    while dates[-1] < end:
        dates.append(dates[-1] + dt.timedelta(days=1))

    # Which CDS dataset to query
    source = "reanalysis-era5-pressure-levels"
    # Shared request parameters
    request = {
        "product_type": "reanalysis",
        "grid": "1.0/1.0",
        "area": "90/-180/0/180",
        "date": ["{:%Y-%m-%d}".format(d) for d in dates],
        "time": ["00:00", "06:00", "12:00", "18:00"],
        "format": "netcdf",
        "variable": [
            "u_component_of_wind",
            "v_component_of_wind",
            "temperature"
        ],
        "pressure_level": [
              "10",   "50",  "100",  "200",  "250",
             "300",  "400",  "500",  "700",  "850",
             "925", "1000"
        ]
    }

    # Name of the output file
    name = "ERA-2016-12-10-to-2016-12-24-uvt.nc"

    # Submit request to the CDSAPI
    c = cdsapi.Client()
    c.retrieve(source, request, name)

