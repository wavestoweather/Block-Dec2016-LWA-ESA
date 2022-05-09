import datetime as dt


start = dt.datetime(2016, 12,  8, 0)
end   = dt.datetime(2016, 12, 18, 0)
valid = dt.datetime(2016, 12, 18, 0)

# FIXME: set an appropriate output path (target parameter)
request = """
retrieve,
    class=od,
    stream=enfo,
    expver=1,
    type=pf,
    grid=1.0/1.0,
    number=1/TO/50,
    area=90/-180/0/180,
    levtype=pl,
    levelist=10/50/100/200/250/300/400/500/700/850/925/1000,
    param=130.128/131/132,
    date={init:%Y-%m-%d},
    time={init:%H}:00:00,
    step={step},
    target="FIXME/ENS-DEC18-EVAL-{init:%Y-%m-%dT%HZ}.grb"
""".format


init = start
while init <= end:
    # Only need data for specified valid time 
    step = int((valid-init).total_seconds() / 60 / 60)
    # Write MARS request
    with open("request-eval-{:%Y-%m-%dT%HZ}".format(init), "w") as f:
        f.write(request(init=init, step=step).lstrip())
    # Next forecast
    init = init + dt.timedelta(hours=12)

