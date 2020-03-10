module Calibration

using Distributions

    function monteCarlo(paramSpace;samples=100000)
        params = Dict{Symbol,Array}()

        keyList = collect(keys(paramSpace))

        for k in keyList
            p = paramSpace[k]
            minV = p[:lower]
            maxV = p[:upper]
            params[k] = collect(rand(Uniform(minV,maxV), samples))
        end

        return params

    end

    function latinHypercube(func,forcing,obs,param_bounds;samples=10000)


    end

# end module
end
