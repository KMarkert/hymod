using CSV, DataFrames, Dates#, Plots
include("src/Hymod.jl")

obsPath = "../example_data/LA011201_obs.csv"
obs = CSV.read(obsPath)
# obs.isodate = Hymod.Utils.parseDates(obs, format="m/d/Y",dateCol=:date)

forcingPath = "../example_data/LA011201_forcings.csv"
forcings = CSV.read(forcingPath)

forcings.pet = Hymod.hargreaves(forcings,tminCol=:tmin,tmaxCol=:tmax,dtCol=:isodate)

paramSpace = Dict{Symbol,Dict}(
    :cmax => Dict{Symbol,Float64}(:lower => 1.0, :upper => 100),
    :bexp => Dict{Symbol,Float64}(:lower => 0.0, :upper => 2.0),
    :alpha => Dict{Symbol,Float64}(:lower => 0.2, :upper => 0.99),
    :ks => Dict{Symbol,Float64}(:lower => 0.01, :upper => 0.5),
    :kq => Dict{Symbol,Float64}(:lower => 0.5, :upper => 1.2)
)

calStart = Date(1986,1,1)
calEnd = Date(1995,12,31)

calForcings = filter(row -> row[:isodate] >= calStart && row[:isodate] <= calEnd, forcings)
obsSel = filter(row -> row[:datetime] >= calStart && row[:datetime] <= calEnd, obs)

nIters = 5

calQ, calPars, loss = Hymod.calibrate(calForcings,obsSel.discharge,paramSpace,nIters)
print(calPars,loss)

# # do some plotting to show the results
# obsSel[:calibrated] = calQ
# obsSel = filter(row -> row[:datetime] >= Date(1987,1,1), obsSel)

# theme(:bright)

# plot(obsSel[:datetime],[obsSel[:discharge] obsSel[:calibrated]],
#     label=["Observed" "Calibrated"],
#     xlabel="Date",
#     xrotation=40,
#     ylabel="Discharge m³/s",
#     dpi=200
# )
# savefig("calibrated_discharge.png")