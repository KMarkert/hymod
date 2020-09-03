import pandas as pd
import datetime
from src.hymod import Hymod

t1 = datetime.datetime.now()
obsPath = "../example_data/LA011201_obs.csv"
obs = pd.read_csv(obsPath)
obs["datetime"] = pd.to_datetime(obs["datetime"])

path = "../example_data/LA011201_forcings.csv"
forcings = pd.read_csv(path)

forcings["datetime"] = pd.to_datetime(forcings["isodate"])
forcings['pet'] = Hymod.hargreaves(forcings,dtCol="datetime")

pars = Hymod.get_random_params()
q = Hymod.simulate(forcings,**pars)
print(f"Processing time: {datetime.datetime.now()-t1}")
