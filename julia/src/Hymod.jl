#
module Hymod

using DataFrames, Dates, Distributions, TOML

include("Utils.jl")
include("Calibration.jl")

    function get_random_params()
        param_bounds = Dict{Symbol,Dict}(
            :cmax => Dict{Symbol,Float64}(:lower => 1.0, :upper => 100),
            :bexp => Dict{Symbol,Float64}(:lower => 0.0, :upper => 2.0),
            :alpha => Dict{Symbol,Float64}(:lower => 0.2, :upper => 0.99),
            :ks => Dict{Symbol,Float64}(:lower => 0.01, :upper => 0.5),
            :kq => Dict{Symbol,Float64}(:lower => 0.5, :upper => 1.2)
        )

        out = Dict()

        for k in keys(param_bounds)
            minV = param_bounds[k][:lower]
            maxV = param_bounds[k][:upper]
            p = rand(Uniform(minV,maxV), 1)
            out[k] = p[1]
        end

        return out
    end

    function _power(x,y)
        x=abs(x) # Needed to capture invalid overflow with netgative values
        return x^y
    end

    function _excess(x_loss,cmax,bexp,prc,pet)
        xn_prev = x_loss
        ct_prev = cmax * (1 - _power((1 - ((bexp + 1) * (xn_prev) / cmax)), (1 / (bexp + 1))))
        # Calculate Effective rainfall 1
        ER1 = max((prc - cmax + ct_prev), 0.0)
        prc = prc- ER1
        dummy = min(((ct_prev + prc) / cmax), 1)
        xn = (cmax / (bexp + 1)) * (1 - _power((1 - dummy), (bexp + 1)))

        # Calculate Effective rainfall 2
        ER2 = max(prc - (xn - xn_prev), 0)

        # Alternative approach
        evap = (1 - (((cmax / (bexp + 1)) - xn) / (cmax / (bexp + 1)))) * pet  # actual ET is linearly related to the soil moisture state
        xn = max(xn - evap, 0)  # update state

        return ER1,ER2,xn
    end

    function linearReservoir(x_slow,inflow,Rs)
        # Linear reservoir routing
        x_slow = (1 - Rs) * x_slow + (1 - Rs) * inflow
        outflow = (Rs / (1 - Rs)) * x_slow
        return x_slow,outflow
    end

    function hargreaves(forcings; tminCol=:tmin,tmaxCol=:tmax,dtCol=:datetime)

        dts = forcings[!,dtCol]
        tmin = forcings[!,tminCol]
        tmax = forcings[!,tmaxCol]
        len = length(tmax)
        Gsc = 367
        lhov = 2.257

        doy = map(dayofyear, dts)

        tavg = map(mean,zip(tmin,tmax))

        eto = zeros(len)

        for (i,t) in enumerate(doy)
            b = 2 * pi * (Int16(t)/365)
            Rav = 1.00011 + 0.034221*cos(b) + 0.00128*sin(b) + 0.000719*cos(2*b) + 0.000077*sin(2*b)
            Ho = ((Gsc * Rav) * 86400)/1e6

            eto[i] = (0.0023 * Ho * (tmax[i]-tmin[i])^0.5 * (tavg[i]+17.8))
        end

        return eto
    end

    function simulate(forcings; precipCol=:precip,petCol=:pet,initFlow=true,cmax=1.0,bexp=0.0,alpha=0.2,ks=0.01,kq=0.5)
        p = forcings[!,precipCol]
        e = forcings[!,petCol]

        lt_to_m = 0.001

        # HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
        x_loss = 0.0
        # Initialize slow tank state
        # value of 0 init flow works ok if calibration data starts with low discharge
        x_slow = initFlow ? (2.3503 / (ks * 22.5)) : 0
        # Initialize state(s) of quick tank(s)
        x_quick = zeros(Float64, 3)
        outflow = zeros(Float64, length(p))
        output = zeros(Float64, length(p))
        # START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS

        for t in 1:length(p)
            Pval = p[t]
            PETval = e[t]
            # Compute excess precipitation and evaporation
            ER1, ER2, x_loss = _excess(x_loss, cmax, bexp, Pval, PETval)
            # Calculate total effective rainfall
            ET = ER1 + ER2
            #  Now partition ER between quick and slow flow reservoirs
            UQ = alpha * ET
            US = (1 - alpha) * ET
            # Route slow flow component with single linear reservoir
            x_slow, QS = linearReservoir(x_slow, US, ks)
            # Route quick flow component with linear reservoirs
            inflow = UQ

            for i = 1:length(x_quick)
                # Linear reservoir
                x_quick[i], outflow = linearReservoir(x_quick[i], inflow, kq)
                inflow = outflow
            end

            # Compute total flow for timestep
            output[t] = ((QS + outflow)/lt_to_m)
        end

        output[output.<0] .= 0

        return output
    end

    function calibrate(forcing, obs, paramSpace, samples; precipCol=:precip, petCol=:pet, saveResults=false)
        # implementation of a monte carlo sampling

        params = Calibration.monteCarlo(paramSpace,samples=samples)
        keyList = collect(keys(params))
        nIter = samples
        losses = zeros(nIter)

        print("Running $nIter iterations...")
        for i in 1:nIter
            vals = [params[k][i] for k in keyList]
            pars = Dict(keyList .=> vals)
            q = simulate(forcing, precipCol=:precip, petCol=:pet; pars...)
            losses[i] = Utils.nse(q[365:end],obs[365:end])
        end

        losses = convert(Vector{Union{Missing,Float64}}, losses)
        losses[map(isnan,losses)] .= missing

        finalLoss, idx = findmax(skipmissing(losses))
        vals = [params[k][idx] for k in keyList]
        finalPars = Dict(keyList .=> vals)
        finalQ = simulate(forcing, precipCol=:precip, petCol=:pet; finalPars...)

        if saveResults
            t = now()
            finalPars[:loss] = finalLoss
            tomlPars = Dict(String(k)=>v for (k,v) in finalPars)
            print(tomlPars)
            open("hymod_calibration_results_$t.toml","w") do io
                TOML.print(io, tomlPars)
            end
        end

        return finalQ, finalPars, finalLoss

    end

# end module
end
