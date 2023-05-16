import numpy
import matplotlib.pyplot as plt
import seaborn
import scipy.stats
import pandas
import astropy.units as units
from astropy import constants
import pickle
import pathlib
from altaipony.ffd import FFD, generate_random_power_law_distribution
g = -1.0
e = 4.39e31


def calculate_log_tau(name):
    grouped = planets.groupby("star_name")
    for idx, group in grouped:
        if idx == name:
            mass_star = group.iloc[0]["star_mass"]
            break
    print(mass_star)
    tau = 0.
    tau_err_min = 0.
    tau_err_max = 0.
    for index, row in logtau.iterrows():
        if row["mass_min"] <= mass_star <= row["mass_max"]:
            print(mass_star)
            tau =  logtau.at[index,"log_tau"]
            print(tau)
            tau_err_min =  logtau.at[index,"log_tau_err_min"]
            tau_err_max =  logtau.at[index,"log_tau_err_max"]

    return tau, tau_err_min, tau_err_max, mass_star


# M dwarfs
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
     "Kepler-94"
 ]
thelistofnamesK2=[
    "K2-240",
    "TRAPPIST-1"
]
thelistofnamesTESS=[
    "TOI-1749"
]

periods = pandas.read_pickle("./data/periods_stars.pkl")
all_flares = pandas.read_pickle("./data/all_flares.pkl")
all_flares = all_flares.query("`ed_ratio`.notna()&`recovery_probability`.notna()&`ed_corr`.notna()")
planets = pandas.read_csv("./data/mdwarf_active_sample.csv",sep=',', index_col=0)
grouped = planets.groupby("star_name")
for idx, group in grouped:
    if idx in periods.index:
        periods.at[idx,"ew_espels_halpha"] = group.iloc[0]["ew_espels_halpha"]

print(periods)
logtau = pandas.read_csv("./data/log_tau.csv",sep=';')
all_rossby = pandas.DataFrame({
        "period_rot": pandas.Series(dtype="float"),
        "freq_from_powerlaw": pandas.Series(dtype="float"),
        "rossby_number": pandas.Series(dtype="float"),
        "ew_espels_halpha":pandas.Series(dtype="float")
    })
for listofstar in [thelistofnamesK2,thelistofnamesTESS,thelistofnamesKepler]:
    for name in listofstar:
            logTau, logTau_err_min, logTau_err_max, mass_star = calculate_log_tau(name)

            print(name)

            rossby = all_flares.query(f"`name` == '{name}'")
            if rossby.empty:
                period_rot = periods.at[name,"period"]
                rossby_val = period_rot/(10 ** logTau)
                all_rossby.at[name,"period_rot"] = period_rot
                all_rossby.at[name,"Mstar"] = mass_star
                all_rossby.at[name,"rossby_number"] = rossby_val
                all_rossby.at[name,"ew_espels_halpha"] = periods.at[name,"ew_espels_halpha"]
                continue
            elif len(rossby) ==1:
                period_rot = periods.at[name,"period"]
                rossby_val = period_rot/(10 ** logTau)
                all_rossby.at[name,"period_rot"] = period_rot
                all_rossby.at[name,"Mstar"] = mass_star
                all_rossby.at[name,"rossby_number"] = rossby_val
                all_rossby.at[name,"ew_espels_halpha"] = periods.at[name,"ew_espels_halpha"]
                continue
            atime = rossby["tstart"].max()-rossby["tstop"].min()

            simple_ffd2 = FFD(f=rossby,ID="name", tot_obs_time=atime)

            simple_ffd2.ed_and_freq(energy_correction=True,
                                   recovery_probability_correction=True)

            simple_ffd2.alpha, simple_ffd2.alpha_err = -g + 1, .1

            simple_ffd2.fit_beta_to_powerlaw()
            rossby = rossby.reset_index()
            print(rossby)
            rate =0.

            if rossby.empty:
                continue
            energy = rossby.at[0,"luminosity quiscent"]
            print(f"luminosity quiscent={energy}")
            ed_cutoff = e/energy
            print(f"ED cutoff={ed_cutoff}, alpha = {simple_ffd2.alpha}, beta = {simple_ffd2.beta}")
            rate = simple_ffd2.beta/ (numpy.abs(simple_ffd2.alpha - 1.))* ed_cutoff ** (1 - simple_ffd2.alpha)
            ed_cutoff_max = rossby["flare_erg_rec"].max()/energy

            period_rot = periods.at[name,"period"]
            rossby_val = period_rot/(10 ** logTau)
            all_rossby.at[name,"period_rot"] = period_rot
            all_rossby.at[name,"Mstar"] = mass_star

            all_rossby.at[name,"freq_from_powerlaw"] = rate
            all_rossby.at[name,"rossby_number"] = rossby_val
            all_rossby.at[name,"ew_espels_halpha"] = periods.at[name,"ew_espels_halpha"]

all_rossby.to_pickle(f"./data/rossby.pkl")

