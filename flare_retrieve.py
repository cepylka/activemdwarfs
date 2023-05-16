from altaipony.lcio import from_mast, from_path
import numpy
import matplotlib.pyplot as plt
import lightkurve as lk
import seaborn
import scipy.stats
import pandas
import astropy.units as units
from astropy import constants
import pickle

response_curve = {"TESS" : "TESS.txt",
                  "Kepler" : "kepler_lowres.txt"}

for key, val in response_curve.items():
    df = pandas.read_csv(f"./data/{val}",
                     delimiter="\s+", skiprows=8)
    df = df.sort_values(by="nm", ascending=True)
    rwav = (df.nm * 10).values * units.angstrom  # convert to angstroms
    rres = (df.response).values
    response_curve[key] = (rwav,rres)

def black_body_spectrum(wav, t):
    """Takes an array of wavelengths and
    a temperature and produces an array
    of fluxes from a thermal spectrum
    Parameters:
    -----------
    wav : Astropy array
        wavenlength array
    t : float
        effective temperature in Kelvin
    """
    t = t * units.K # set unit to Kelvin

    return (( (2 * numpy.pi * constants.h * constants.c**2) / (wav**5) / (numpy.exp( (constants.h * constants.c) / (wav * constants.k_B * t) ) - 1))
            .to("erg*s**(-1)*cm**(-3)")) #simplify the units

def calculate_specific_flare_flux(mission, flaret=1e4):
    """ Get specific Kepler/TESS flux
    Parameters:
    -----------
    mission : string
        TESS or Kepler
    flaret : float
        black body temperature
    Return:
    -------
    specific Kepler/TESS flux in units erg*cm**(-2)*s**(-1)
    """

    try:
        # Read in response curve:
        rwav, rres = response_curve[mission]
    except KeyError:
        raise KeyError("Mission can be either Kepler or TESS.")
    # create an array to upsample the filter curve to
    w = numpy.arange(3000,13001) * units.angstrom

    # interpolate thermal spectrum onto response
    # curve wavelength array, then sum up
    # flux times response curve:

    # Generate a thermal spectrum at the Teff
    # over an array of wavelength w:
    thermf = black_body_spectrum(w, flaret)

    # Interpolate response from rwav to w:
    rress = numpy.interp(w,rwav,rres)#, left=0, right=1)

    # Integrating the flux of the thermal
    # spectrum times the response curve over wavelength:
    return numpy.trapz(thermf * rress, x=w).to("erg*cm**(-2)*s**(-1)")


planets = pandas.read_csv("./data/mdwarf_active_sample.csv",sep=',', index_col=0)

i_units = units.Quantity(1, unit="erg cm-2 s-1")


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
    # "Kepler-114",
    # "Kepler-138", #use cadence="fast"
    # "Kepler-32",
    # "Kepler-446",
    # "Kepler-52",
    # "Kepler-83",
    # "Kepler-94",
    "TOI-1749"
]

for listofstar in [thelistofnamesK2,thelistofnamesTESS,thelistofnamesKepler]:

    if listofstar == thelistofnamesKepler:
        lmt = 10
        for name in listofstar:
            print(name)
            kepler_flux = calculate_specific_flare_flux("Kepler", flaret=planets.at[f"{name} b","star_teff"])
            # quiscent star luminosity erg/s
            luminosity = 4 * kepler_flux * numpy.pi * (planets.at[f"{name} b","star_radius"]*constants.R_sun.to(units.cm)) ** 2

            all_flares = pandas.DataFrame()
            quarters = pandas.DataFrame({
                    "t_min": pandas.Series(dtype="float"),
                    "t_min": pandas.Series(dtype="float"),
                })
            while True:
                try:
                    flc = from_mast(
                        name,
                        mode="LC",
                        cadence="short",
                        mission="Kepler"
                    )
                except Exception as ex:
                    print(f"Exception: {type(ex)}")
                    lmt -= 1
                    if lmt <= 0:
                        break # end of download's tries
                    else:
                        continue
                break # end of while true cycle


            for j in range(len(flc)):
                t_min = flc[j].time.min()
                t_max = flc[j].time.max()
                quarters.loc[j,"t_min"] = t_min.value
                quarters.loc[j,"t_max"] = t_max.value
                quarters["name"] = name
                if flc[j] == None:
                    continue

                flcd = flc[j].detrend("savgol") # detrend curve
                flcd = flcd.find_flares(N1=3, N2=1, N3=3,minsep=3) # find flares
                flcdpanda = flcd.flares
                if not flcdpanda.empty:
                    print(flcdpanda)
                    flcd, fakeflc = flcd.sample_flare_recovery(inject_before_detrending=True, mode="savgol",
                                                  iterations=20, fakefreq=2, ampl=[1e-4, 0.5],
                                                   dur=[.001/6., 0.1/6.]
                                                  ) # injection of simulated flares
                    print("The total number of injected flares is {}.".format(flcd.fake_flares.shape[0]))
                    flcc = flcd.characterize_flares(ampl_bins=10, dur_bins=10)
                    flcdpanda = flcc.flares
                    flcdpanda["luminosity quiscent"] = luminosity

                    flcdpanda["flare_erg"] = flcdpanda["ed_rec"]*luminosity
                    flcdpanda["flare_erg_rec"] = flcdpanda["ed_corr"]*luminosity
                    print(flcdpanda)
                    all_flares = pandas.concat([all_flares,flcdpanda],ignore_index=True)
                else:
                    print('DataFrame is empty! No flares')

            all_flares.to_pickle(f"./data/flare_pkl/{name}_flares.pkl") # table of flare parameters
            quarters.to_pickle(f"./data/quarters_pkl/{name}_quarters.pkl") # table of time stamps of quarters with flares

    elif listofstar == thelistofnamesTESS:
        lmt = 10
        for name in listofstar:
            print(name)
            tess_flux =  calculate_specific_flare_flux("TESS", flaret=planets.at[f"{name} b","star_teff"])
            # quiscent star luminosity erg/s
            luminosity_tess = 4 * tess_flux * numpy.pi * (planets.at[f"{name} b","star_radius"]*constants.R_sun.to(units.cm)) ** 2
            luminosity = luminosity_tess/0.72
            all_flares = pandas.DataFrame()
            quarters = pandas.DataFrame({},index=None)


            while True:
                try:
                    flc = from_mast(
                        name,
                        mode="LC",
                        cadence="short",#"fast", # if to change for 20 s cadence
                        mission="TESS",
                        author="SPOC"
                    )
                except Exception as ex:
                    print(f"Exception: {type(ex)}")
                    lmt -= 1
                    if lmt <= 0:
                        break
                    else:
                        continue
                break
            for j in range(len(flc)):
                t_min = flc[j].time.min()
                t_max = flc[j].time.max()

                quarters.at[j,"t_min"] = t_min.value
                quarters.at[j,"t_max"] = t_max.value
                quarters["name"] = name
                if flc[j] == None:
                    continue

                flcd = flc[j].detrend("savgol") # detrend curve
                flcd = flcd.find_flares(N1=3, N2=1, N3=3,minsep=3) # find flares
                flcdpanda = flcd.flares
                if not flcdpanda.empty:
                    flcd, fakeflc = flcd.sample_flare_recovery(inject_before_detrending=True, mode="savgol",
                                                  iterations=30, fakefreq=1, ampl=[1e-4, 0.5],
                                                   dur=[.001/6., 0.1/6.]
                                                  ) # injection of simulated flares
                    print("The total number of injected flares is {}.".format(flcd.fake_flares.shape[0]))
                    flcc = flcd.characterize_flares(ampl_bins=10, dur_bins=10)
                    flcdpanda = flcc.flares
                    flcdpanda["luminosity quiscent"] = luminosity
                    flcdpanda["flare_erg"] = flcdpanda["ed_rec"]*luminosity
                    flcdpanda["flare_erg_rec"] = flcdpanda["ed_corr"]*luminosity
                    print(flcdpanda)
                    all_flares = pandas.concat([all_flares,flcdpanda],ignore_index=True)
                else:
                    print('DataFrame is empty! No flares')

            all_flares.to_pickle(f"./data/flare_pkl/{name}_flares.pkl") # table of flare parameters
            quarters.to_pickle(f"./data/quarters_pkl/{name}_quarters.pkl") # table of time stamps of quarters with flares

    else:
        for name in listofstar:
            print(name)
            kepler_flux = calculate_specific_flare_flux("Kepler", flaret=planets.at[f"{name} b","star_teff"])
            # quiscent star luminosity erg/s
            luminosity = 4 * kepler_flux * numpy.pi * (planets.at[f"{name} b","star_radius"]*constants.R_sun.to(units.cm)) ** 2

            all_flares = pandas.DataFrame()
            quarters = pandas.DataFrame({
                    "t_min": pandas.Series(dtype="float"),
                    "t_min": pandas.Series(dtype="float"),
                })
            lmt = 10
            if name == "K2-240":
                while True:
                    try:
                        flc = from_mast(name,
                            mode="LC",
                            cadence="long",
                            mission="K2",
                        ) # only LC curves for K2-240
                    except Exception as ex:
                        print(f"Exception: {type(ex)}")
                        lmt -= 1
                        if lmt <= 0:
                            break
                        else:
                            continue
                    break

                t_min = flc.time.min()
                t_max = flc.time.max()
                print(t_min,t_max)
                quarters.loc[0,"t_min"] = t_min.value
                quarters.loc[0,"t_max"] = t_max.value
                quarters["name"] = name
                flcd = flc.detrend("savgol")
                flcd = flcd.find_flares(N1=3, N2=1, N3=3,minsep=3)
                flcdpanda = flcd.flares
                flcdpanda = flcd.flares
                if not flcdpanda.empty:
                    print(flcdpanda)
                    flcd, fakeflc = flcd.sample_flare_recovery(inject_before_detrending=True, mode="savgol",
                                                  iterations=30, fakefreq=3, ampl=[1e-4, 0.5],
                                                   dur=[.001/6., 0.1/6.]
                                                  )
                    print("The total number of injected flares is {}.".format(flcd.fake_flares.shape[0]))
                    flcc = flcd.characterize_flares(ampl_bins=10, dur_bins=10)
                    flcdpanda = flcc.flares
                    flcdpanda["luminosity quiscent"] = luminosity
                    flcdpanda["flare_erg"] = flcdpanda["ed_rec"]*luminosity
                    flcdpanda["flare_erg_rec"] = flcdpanda["ed_corr"]*luminosity
                    print(flcdpanda)
                    all_flares = pandas.concat([all_flares,flcdpanda],ignore_index=True)
                else:
                    print('DataFrame is empty! No flares')
            elif name == "TRAPPIST-1":
                lmt = 10

                while True:
                    try:
                        flc = from_mast(name,
                            mode="LC",
                            cadence="short",#"long"
                            mission="K2",
                        )
                    except Exception as ex:
                        print(f"Exception: {type(ex)}")
                        lmt -= 1
                        if lmt <= 0:
                            break
                        else:
                            continue
                    break

                t_min = flc[0].time.min()
                t_max = flc[0].time.max()
                quarters.at[0,"t_min"] = t_min.value
                quarters.at[0,"t_max"] = t_max.value
                quarters["name"] = name
                flcd = flc[0].detrend("savgol")
                flcd = flcd.find_flares(N1=3, N2=1, N3=3,minsep=3)
                flcdpanda = flcd.flares
                if not flcdpanda.empty:
                    print(flcdpanda)
                    flcd, fakeflc = flcd.sample_flare_recovery(inject_before_detrending=True, mode="savgol",
                                                  iterations=30, fakefreq=1, ampl=[1e-4, 0.5],
                                                   dur=[.001/6., 0.1/6.]
                                                  )
                    print("The total number of injected flares is {}.".format(flcd.fake_flares.shape[0]))
                    flcc = flcd.characterize_flares(ampl_bins=10, dur_bins=10)
                    flcdpanda = flcc.flares
                    flcdpanda["luminosity quiscent"] = luminosity
                    flcdpanda["flare_erg"] = flcdpanda["ed_rec"]*luminosity
                    flcdpanda["flare_erg_rec"] = flcdpanda["ed_corr"]*luminosity
                    print(flcdpanda)
                    all_flares = pandas.concat([all_flares,flcdpanda],ignore_index=True)
                else:
                    print('DataFrame is empty! No flares')

            all_flares.to_pickle(f"./data/flare_pkl/{name}_flares.pkl")
            quarters.to_pickle(f"./data/quarters_pkl/{name}_quarters.pkl")

