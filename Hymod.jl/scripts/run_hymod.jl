using CSV, DataFrames, Dates, Gadfly
include("../src/Hymod.jl")

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
calEnd = Date(1990,12,31)

calForcings = filter(row -> row[:isodate] >= calStart && row[:isodate] <= calEnd, forcings)
obsSel = filter(row -> row[:datetime] >= calStart && row[:datetime] <= calEnd, obs)

calQ, calPars, loss = Hymod.calibrate(calForcings,obsSel.discharge,paramSpace,1000000)
print(calPars,loss)
#
obsSel[:calibrated] = calQ
obsSel = filter(row -> row[:datetime] >= Date(1987,1,1), obsSel)

p = plot(obsSel,layer(x=:datetime,y=:calibrated,Geom.line,Theme(default_color=color("orange"))),
                layer(x=:datetime,y=:discharge, Geom.line,Theme(default_color=color("purple"))))

img = SVG("hymod_dischage.svg", 20cm, 16cm)
draw(img, p)
