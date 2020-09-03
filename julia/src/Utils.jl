module Utils

using DataFrames, Dates, Statistics

    function parseDates(df; format="y-m-d",dateCol=:date)
        dfmat = DateFormat(format)
        tCol = df[!,dateCol]
        len = length(tCol)
        dts = Array{Date}(undef,len)
        for i in 1:len
            dts[i] = Date(tCol[i],dfmat)
        end
        return dts
    end

    function renameCols(df,newNames)
        @assert length(names(df)) == length(newNames) "List of new names does not equal to columns"
        names!(df,[Symbol(i) for i in newNames])
        return df
    end

    function nse(sim,obs)
        numerator = sum((obs-sim).^2)
        denominator = sum((obs.-mean(obs)).^2)
        return 1 - (numerator/denominator)
    end

end
