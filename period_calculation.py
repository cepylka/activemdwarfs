import numpy
import matplotlib.pyplot as plt
import lightkurve as lk
import pandas
import seaborn
import scipy.stats



thelistofnamesKepler=[
    "Kepler-114",
    "Kepler-138",
    "Kepler-231",
    "Kepler-26",
    "Kepler-28",
    "Kepler-32",
    "Kepler-446",
    "Kepler-49",
    "Kepler-52",
    "Kepler-83",
    "Kepler-94",
 ]
thelistofnamesK2=[
    "K2-240",
    "TRAPPIST-1",
]
thelistofnamesTESS=[
    "TOI-1749",
    # "Kepler-114",
    # "Kepler-138",
    # "Kepler-32",
    # "Kepler-446",
    # "Kepler-52",
    # "Kepler-83",
    # "Kepler-94",
]
periods_stars = pandas.DataFrame()
all_stars = []
for star in thelistofnamesKepler:
    all_stars.append(star)
for star1 in thelistofnamesTESS:
    all_stars.append(star1)
for star2 in thelistofnamesK2:
    all_stars.append(star2)
print(list(all_stars))
periods_stars["period"] = numpy.array(numpy.NaN, dtype=float)
periods_stars["uid"] = all_stars
periods_stars["uid"] = periods_stars["uid"].astype('str')

print(periods_stars["uid"])
periods_stars.set_index('uid',inplace = True)

for listofstar in [thelistofnamesKepler,thelistofnamesTESS,thelistofnamesK2]:
    if listofstar == thelistofnamesKepler:
        for name in listofstar:
            period_mean = 0.
            print(name)
            if name == "Kepler-114":
                search_result = lk.search_lightcurve(
                "kplr010925104",
                author='Kepler',
                cadence="long") # the standart search additionaly returns one other star, to avoid it we state the name as in MAST cataloque
            else:
                search_result = lk.search_lightcurve(
                    name,
                    author='Kepler',
                    cadence="long")

            lc = search_result.download_all().stitch()
            no_nan_lc = lc.remove_nans()
            clipped_lc = no_nan_lc.remove_outliers(sigma=3)
            pg = clipped_lc.to_periodogram(maximum_period=100)
            period = pg.period_at_max_power
            period = period.value
            print(period)
            periods_stars.at[name,"period"] = period
    elif listofstar == thelistofnamesTESS:
        for name in listofstar:
            period_mean = 0.
            print(name)
            search_result = lk.search_lightcurve(
            name,
            author='TESS-SPOC',
            cadence="long")
            lc = search_result.download_all().stitch()
            no_nan_lc = lc.remove_nans()
            clipped_lc = no_nan_lc.remove_outliers(sigma=3)
            pg = clipped_lc.to_periodogram(maximum_period=100)
            period = pg.period_at_max_power
            period = period.value
            print(period)
            periods_stars.at[name,"period"] = period

    else:
        for name in listofstar:
            period_mean = 0.
            print(name)
            search_result = lk.search_lightcurve(
            name,
            author="K2",
            cadence="long")

            lc = search_result.download_all().stitch()
            no_nan_lc = lc.remove_nans()
            clipped_lc = no_nan_lc.remove_outliers(sigma=3)
            pg = clipped_lc.to_periodogram(maximum_period=100)
            period = pg.period_at_max_power
            period = period.value
            print(period)
            periods_stars.at[name,"period"] = period

print(periods_stars)
periods_stars.to_pickle("./data/periods_stars.pkl")