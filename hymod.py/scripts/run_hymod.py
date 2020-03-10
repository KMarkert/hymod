import pandas as pd
from src.hymod import Hymod

obsPath = "../data/LA011201_obs.csv"
obs = pd.read_csv(obsPath)
obs["datetime"] = pd.to_datetime(obs["date"])

path = "../data/LA011201_forcings.csv"
forcings = pd.read_csv(path)

forcings["datetime"] = pd.to_datetime(forcings["isodate"])
forcings['pet'] = Hymod.hargreaves(forcings,dtCol="datetime")

paramSpace = dict(
    cmax = dict(lower = 1.0, upper = 100, samples = 10),
    bexp = dict(lower = 0.0, upper = 2.0, samples = 10),
    alpha = dict(lower = 0.2, upper = 0.99, samples = 10),
    ks = dict(lower = 0.01, upper = 0.5, samples = 10),
    kq = dict(lower = 0.5, upper = 1.2, samples = 10)
)

calStart = "1981-01-01"
calEnd = "1984-12-31"

forcingMask = (forcings["datetime"]>= calStart) & (forcings["datetime"]<=calEnd)
calForcings = forcings.loc[forcingMask]

obsMask = (obs["datetime"]>= calStart) & (obs["datetime"]<=calEnd)
obsQ = obs.loc[obsMask]["discharge"]

calQ, calPars, loss = Hymod.calibrate(calForcings,obsQ,paramSpace)

print(calPars, loss)
