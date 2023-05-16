import pandas
import pathlib

pkls = sorted(pathlib.Path("./data/flare_pkl/").rglob("*.pkl"))
all_flares = pandas.DataFrame()
j=0
for pkl in pkls:
    p = pkl
    table = pandas.DataFrame()
    table = pandas.read_pickle(pkl)

    table["name"] = p.stem.split("_")[0]
    all_flares = pandas.concat([all_flares,table])
    j+=1
all_flares = all_flares.reset_index()
print(all_flares)
all_flares = all_flares.drop(["index"], axis=1)

print(all_flares)
all_flares.to_pickle("./data/all_flares.pkl")