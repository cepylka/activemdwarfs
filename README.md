# Active M dwarf with exoplanetary systems

## Preparing the sample

For retrieving information from working databases initially we need a list of potentially interesting systems. We obtained a sample from [Mamonova, in prep.](#Mamonova). For this data we retrieved a list of potential targets, stars of M spectral class for the next step.

## Retrieving stellar parameters from GAIA Archive

The Python module `uio.task` from [uio-exoplanet-group](https://github.com/retifrav/uio-exoplanet-group) (see the documentation [here](https://uio.decovar.dev/uio/tasks/reconfirming_stellar_parameters.html)) is used for retrieving certain stellar parameters from GAIA DR3 release [Fouesneau et al., 2022](#Fouesneau).

``` py
from uio.tasks import reconfirming_stellar_parameters
```

In this case, the module works with the preexisting file `mdwarf_active_sample.pkl` and, from the table `gaiadr3.astrophysical_parameters`, it retrieves, for example, `ew_espels_halpha` and `ew_espels_halpha_model` parameters. You can find this code in the Systems.ipynb. Please use pandas=>2.0.1 for this specific notebook.

``` py
tbl = reconfirming_stellar_parameters.lookForParametersInGaia(
    "./data/mdwarf_active_sample.pkl",
    "gaiadr3.astrophysical_parameters",
    [
        "ew_espels_halpha",
        "ew_espels_halpha_model"
    ],
    "dr3"
)
```
For an explanation of these parameters, we refer to [Creevey et al., 2022](#Creevey) and [Fouesneau et al., 2022](#Fouesneau). Then the module writes the retrieved data in the same table, which can be saved in .pkl file.

After the selection of active M dwarf from the list of planets, performed by analysing found `ew_espels_halpha` parameter, we prepare a new file `mdwarf_active_sample.pkl` comprising parameters of systems in our sample.

## Sample star periods

In order to analyse the sample's stars, you would need to install the Python module `lightkurve` (https://docs.lightkurve.org/about/citing.html).

``` sh
$ pip install lightkurve
```

We use `period_calculation.py` for determination of rotation periods of stars in the sample. You could run this script in the terminal.

``` sh
$ python ./period_calculation.py
```

It will write the file `periods_stars.pkl` to your `./data` folder in the current directory.

For finding periods of Kepler stars using TESS data, you should edit the script: comment all Kepler names in `thelistofnamesKepler` and uncomment them in `thelistofnamesTESS`. Do not forget to change the folder where you want to save new .pkl.

Visualising of period retrieving was also implemented in the jupyter notebook `./notebooks/Mdwarfs_periods.ipynb`

In order to analyse the flares the sample's stars exhibited during the periods of observation, you would need to install the Python module `altaipony` ([Ilin et al. (2021)](#A) and [Davenport (2016)](#B)). Make sure that your system meets all requirements we stated in the file `requirments_altaipony.txt` , otherwise it could not work properly.

``` sh
$ pip install altaipony
```
We use `flare_retrieve.py` for flare finding in the sample's stars. You could run this script in the terminal.

``` sh
$ python ./flare_retrieve.py
```

The script works with our sample stars, and for each of them detrends individual light curves, finds flares, and performs the injection-recovery procedure to obtain corrected parameters and the probability of recovering for each flare event.
The results will be written in your data folder for each star separately: `./data/flare_pkl/{name}_flares.pkl` for flares data (`name` is a star identifier). The time data from quaters, during wich the flares were observed, will be stored in files `./data/quaters_pkl/{name}_quarters.pkl`. The flare-finding and injection-recovering procedures are time-consuming and it would take a couple of hours to run this script.

For finding flares in Kepler stars using TESS data, you should edit the script: comment all Kepler names in `thelistofnamesKepler` and uncomment them in `thelistofnamesTESS`. Do not forget to change the folder where you want to save new .pkl.

``` py
all_flares.to_pickle(f"./data/flare_pkl/your_folder/{name}_flares.pkl")

```

We found flares in TESS data for Kepler-138 only, both with `cadence="short"` and `cadence="fast"` modes. Flare energies would be converted to Kepler energies by default, so if you want to combine all found flares in one Pandas data frame, simply use `pandas.concat`:

``` py
tbl = pandas.read_pickle("./data/flare_pkl/{name}_flares.pkl")
tbl1 = pandas.read_pickle("./data/flare_pkl/your_folder/{name}_flares.pkl")

tbl_new = pandas.concat([tbl, tbl1])
tbl_new.reset_index()
tbl2 = pandas.read_pickle("./data/quarters_pkl/{name}_flares.pkl")
tbl3 = pandas.read_pickle("./data/quarters_pkl/your_folder/{name}_flares.pkl")

tbl_new1 = pandas.concat([tbl2, tbl3])
tbl_new1.reset_index()
```

For subsequent analyses you would need all flares from all stars combined in one pandas data frame using the script below:

``` sh
$ python ./combine_flares.py
```

This script will write `all_flares.pkl` in your data folder. Now you have the data for analyses, and you can use Jupyter notebooks in `notebooks` folder for that.

The script `Rossby_calculation.py` works with the preexisting file `./data/mdwarf_active_sample.pkl`, obtained star periods in file `./data/periods_stars.pkl`, and file `./data/all_flares.pkl` containing all flares data. You would also need file `./data/log_tau.csv`, which is log10 convective turnover times for stars with different masses, reported by [Wright et al. 2018](#Wright). The script calculates star-wise: Rossby numbers, flare frequency at energy cut-off log10E=31.5,

``` sh
$ python ./Rossby_calculation.py
```

You can use the Jupyther notebook `Comparing_samples.ipynb` for comparative analyses. The other sample is from [Medina et al. 2020](#Medina).

## References

- <a name="A"></a> Ekaterina Ilin, Sarah J. Schmidt, Katja Poppenhäger, James R. A. Davenport, Martti H. Kristiansen, Mark Omohundro (2021). "Flares in Open Clusters with K2. II. Pleiades, Hyades, Praesepe, Ruprecht 147, and M67" Astronomy & Astrophysics, Volume 645, id.A42, 25 pp. https://doi.org/10.1051/0004-6361/202039198
- <a name="B"></a> James R. A. Davenport "The Kepler Catalog of Stellar Flares" The Astrophysical Journal, Volume 829, Issue 1, article id. 23, 12 pp. (2016). https://doi.org/10.3847/0004-637X/829/1/23
- <a name="Andrae"></a>Andrae, R. et al. (2022) ‘Gaia Data Release 3: Analysis of the Gaia BP/RP spectra using the General Stellar Parameterizer from Photometry’, Astronomy & Astrophysics [Preprint]. Available at: https://doi.org/10.1051/0004-6361/202243462.
- <a name="Creevey"></a>Creevey, O.L. et al. (2022) ‘Gaia Data Release 3: Astrophysical parameters inference system (Apsis) I -- methods and content overview’. arXiv. Available at: http://arxiv.org/abs/2206.05864 (Accessed: 13 February 2023).
- <a name="Fouesneau"></a>Fouesneau, M. et al. (2022) ‘Gaia Data Release 3. Apsis II: Stellar parameters’, Astronomy & Astrophysics [Preprint]. Available at: https://doi.org/10.1051/0004-6361/202243919.
- <a name="Mamonova"></a>Mamonova, E. et al. (no date) ‘Patterns in the sky. Limited similarity in exoplanetary systems.’
- <a name="Medina"></a>Medina, A.A., Winters, J.G., Irwin, J.M. and Charbonneau, D., 2020. Flare Rates, Rotation Periods, and Spectroscopic Activity Indicators of a Volume-complete Sample of Mid-to Late-M Dwarfs within 15 pc. The Astrophysical Journal, 905(2), p.107. Available at: https://doi.org/10.3847/1538-4357/abc686
- <a name="Wright"></a> Wright, N.J., Newton, E.R., Williams, P.K., Drake, J.J. and Yadav, R.K., 2018. The stellar rotation–activity relationship in fully convective M dwarfs. Monthly Notices of the Royal Astronomical Society, 479(2), pp.2351-2360. Available at: https://doi.org/10.1093/mnras/sty1670.
- <a name="Medina"></a>Medina, A.A., Winters, J.G., Irwin, J.M. and Charbonneau, D., 2020. Flare Rates, Rotation Periods, and Spectroscopic Activity Indicators of a Volume-complete Sample of Mid-to Late-M Dwarfs within 15 pc. The Astrophysical Journal, 905(2), p.107. Available at: https://doi.org/10.3847/1538-4357/abc686
