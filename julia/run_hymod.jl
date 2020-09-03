using CSV, DataFrames, Dates, Plots
include("src/Hymod.jl")

@time begin
    obsPath = "../example_data/LA011201_obs.csv"
    obs = CSV.read(obsPath)
    # obs.isodate = Hymod.Utils.parseDates(obs, format="m/d/Y",dateCol=:date)

    forcingPath = "../example_data/LA011201_forcings.csv"
    forcings = CSV.read(forcingPath)

    forcings.pet = Hymod.hargreaves(forcings,tminCol=:tmin,tmaxCol=:tmax,dtCol=:isodate)

    pars = Hymod.get_random_params()
    q = Hymod.simulate(forcings,precipCol=:precip, petCol=:pet; pars...)
end
