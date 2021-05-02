using NetCDF, Statistics
using Plots, StatsPlots, ColorSchemes

##### Namelists
myPath = "/home/xinhou/radPV/"
vars = ["tas", "rsds", "rsdscs", "clt", "prw", "pv"]
ssps = ["historical", "ssp126", "ssp245", "ssp370", "ssp585"];

##### functions
function rmcp(path)
    cd(path)
    rm(".ipynb_checkpoints", recursive=true, force=true)
    rm(".directory", recursive=true, force=true)
end

avg(arr; n=3) = dropdims(mean(arr; dims=n); dims=n) 

function get_percentile(arr, pct) # array, percent
    r = zeros(size(arr,1))
    for i = 1:size(arr,1)
        v = arr[i,:]
        ind = floor(Int64, pct*length(v))
        narr = sort(v);
        r[i] = narr[ind];
    end
    return r
end

function rib(arr, md) # array, median, ribbon between the percentiles
    mx = get_percentile(arr, 0.8)
    mn = get_percentile(arr, 0.2)
    return (md.-mn, mx.-md)
end

pathPv = joinpath(myPath, "pv", ssps[1], "day/EU/regr")
rmcp(pathPv)
fnGr = joinpath(pathPv, readdir(pathPv)[end])
area = ncread(fnGr, "cell_area") ./ 1e10;

function getArr(var, s, m) # get array by variable, ssp, model
    path = joinpath(myPath, var, ssps[s], "day/EU/tilt_lat/regr")
    rmcp(path)
    fn = joinpath(path, readdir(path)[m])
    arr = ncread(fn, var)
    return arr
end

###### day indices of various models and years
days = 365 *ones(Int64, 20,)
daysLp = copy(days)
# historical leap years
daysLp[[2,6,10,14,18]] .= 366 # [1996, 2000, 2004, 2008, 2012]
daysLpIdxEnd = [sum(daysLp[1:y]) for y in 1:20]
daysLpIdxHd = daysLpIdxEnd .- daysLp .+1
# ssp
daysLp1 = copy(days)
daysLp1[[4,8,12,16]] .= 366 # [2084. 2088, 2092, 2096]
daysLp1IdxEnd = [sum(daysLp1[1:y]) for y in 1:20]
daysLp1IdxHd = daysLp1IdxEnd .- daysLp1 .+1;

function getPvMatrix(arr) # yield 20-yr series of 1 grid in shape of days x 20 (yr)
    len = length(arr)
    if len % 20 != 0 # 7305, 7304
        if len == 7305
            IdxHd, IdxEnd = daysLpIdxHd, daysLpIdxEnd
        else
            IdxHd, IdxEnd = daysLp1IdxHd, daysLp1IdxEnd
        end
        ar1 = arr[IdxHd[1]:IdxEnd[1]]
        for y in 2:20
            ar1dNew = arr[IdxHd[y]:IdxEnd[y]]
            if length(ar1dNew) == 366
                deleteat!(ar1dNew,60) # del Feb.29
            end
            ar1 = hcat(ar1, ar1dNew)
        end
    else # 7300, 7160
        ar1 = reshape(arr, (Int(len/20), 20))
    end    
    return ar1
end

##### Geo-maps
using PyCall, PyPlot

plt = pyimport("matplotlib.pyplot")
clrs = pyimport("matplotlib.colors")
ccrs = pyimport("cartopy.crs")
cput = pyimport("cartopy.util")
cfeature = pyimport("cartopy.feature")

proj1 = ccrs.PlateCarree();
proj2 = ccrs.Robinson();
proj3 = ccrs.EuroPP();

function agreeDots(ch) # calculate if less than 3/4 of models agree on the sign of change
    lolim = round(Int, size(ch, 3) * 1/4) # lower limit
    uplim = round(Int, size(ch, 3) * 3/4) # upper limit

    ar = lolim .< dropdims(sum(ch .> 0, dims=3); dims=3) .< uplim
    arr = replace(ar, 0 => NaN)
    return arr'
end