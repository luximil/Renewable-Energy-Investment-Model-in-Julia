# This is Julia program written by
# - Francisco Jose Manjon Cabeza Garcia,
# - Alessandro James Carenini,
# - Francesca Cioccarelli
# - and Gaspard Robert
# to implement the investment optimisation model described in our term paper for the energy modelling
# course at TU Berlin in the summer semester 2022. The code is mostly based on the code provided by
# the course lecturers in exercise session 5, in the Julia file called 05-invest01.jl.
#
# This optimisation model returns the optimal investments in dispatchable, storage and indispatchable
# electricity production facilities, as well as the optimal irrigation strategy for a farm to reduce the
# total cost of investments plus electricity production costs.

using JuMP # building models.
using JuMP.Containers # use JuMP containers for parameters.
using GLPK # solver for the JuMP model.
using CSV # readin of CSV files.
using DataFrames # data tables.
using CairoMakie # results plotting.

# -----------------------------------------------------------------------------------------------------------
# User parameters                                                                                           |
# -----------------------------------------------------------------------------------------------------------
# Define here the path to the folder containing the model parameters in csv files.
data_path = joinpath("Data")


# -----------------------------------------------------------------------------------------------------------
# Data imports                                                                                              |
# -----------------------------------------------------------------------------------------------------------

# DataFrame with data for temperature, precipitation, humidity, solar energy, solar radiation, etc.
# for each time step. This data file can be obtained at https://www.visualcrossing.com/weather-history/
# The model needs the data from the following columns (create a Weather_data.csv file with those columns and
# column names, if you are not using the data from VisualCrossing):
# - temp:               (float) temperature in degrees Celsius.
# - precip:             (float) precipitation in mm.
# - windspeed:          (float) wind speed in km/h.
# - solarradiation:     (float) solar radiation.
# - solarenergy:        (float) solar energy.
df_weather_time = CSV.read(joinpath(data_path, "Weather_data.csv"), DataFrame)

# DataFrame with parameters for supply side technologies: dispatchable, storage and indispatchable facilities.
# The data frame contains the following columns:
# - technology:         (string) name of the technology. PV and wind facility's names must start with "pv" or
#                       "wind", so that the program can identify them as such and use the correct availability
#                       parameter in the computations.
# - dispatchable:       (bool) if the facility is dispatchable or not.
# - renewable:          (bool) if the facility produces renewable energy or not.
# - variable_cost:      (float) variable cost of electricity production per unit of capacity.
# - fixed_cost:         (float) fixed cost of the facility per unit of capacity.
# - eff:                (float) efficiency factor.
# - inv_power:          (float) investment cost per unit of capacity for electricity producing facilities. Set
#                       this parameter to 0 for storage facilites.
# - inv_storage:        (float) investment cost per unit of capacity for storage facilities. Set this parameter
#                       to 0 for non-storage facilites.
# - lifetime:           (float) lifetime of the facility.
# - emission_factor:    (float) emission factor of the facility.
# - fuel_cost:          (float) fuel cost per unit of produced electricity.
# - surface_covered:    (float) surface covered by the facility per unit of capacity.
# - max_surface:        (float) maximum surface constraint of the facility.
# - g_max:              (float) if set to a number strictly greater than 0, the facility is considered to be
#                       installed, and its capacity is not optimised by the model, i.e., is assumed to be fixed.
df_supply_tech = CSV.read(joinpath(data_path, "Supply_technology.csv"), DataFrame)

# DataFrame with the parameters of the water side technologies: water pumps and water tanks. The DataFrame
# contains the following columns:
# - technology:         (string) name of the technology.
# - storage:            (bool) if the facility is a water tank or not.
# - energy_cost:        (float) cost of pumping one unit of water to the field in units of electricity.
# - max_capacity:       (float) maximum capacity if the facility is a water tank.
df_water_tech = CSV.read(joinpath(data_path, "Water_technology.csv"), DataFrame)

# DataFrame with the hourly electricity demand of each of the other demand sources as columns.
df_other_demand = CSV.read(joinpath(data_path, "Other_demand_data.csv"), DataFrame)

# DataFrame with other parameters. The DataFrame contains the following rows:
# - CO2_price:          (float) CO2 price.
# - Lost_load_cost:     (float) lost load cost.
# - Interest_rate:      (float) interest rate.
# - Daily_effective_water (float) daily effective water need of the field.
# - evtr_alpha:         (float) alpha parameter for the evapotranspiration computation.
# - altitude:           (float) altitude parameter for the evapotranspiration computation.
df_other_parameters = CSV.read(joinpath(data_path, "Other_parameters.csv"), DataFrame)

# Save the CO2 price to a variable.
co2_price = float.(df_other_parameters[df_other_parameters.Parameter .== "CO2_price", "Value"])[1]
# Save the lost load cost to a variable.
lost_load_cost = float.(df_other_parameters[df_other_parameters.Parameter .== "Lost_load_cost", "Value"])[1]
# Save the interest rate to a variable.
interest_rate = float.(df_other_parameters[df_other_parameters.Parameter .== "Interest_rate", "Value"])[1]
# Save the daily effective irrigation needed by the field to a variable.
DEW = float.(df_other_parameters[df_other_parameters.Parameter .== "Daily_effective_water", "Value"])[1]
# Save the alpha parameter for the evapotranspiration computation.
alpha = float.(df_other_parameters[df_other_parameters.Parameter .== "evtr_alpha", "Value"])[1]
# Save the altitude parameter for the evapotranspiration computation.
altitude = float.(df_other_parameters[df_other_parameters.Parameter .== "altitude", "Value"])[1]

# -----------------------------------------------------------------------------------------------------------
# Model parameter definitions                                                                               |
# -----------------------------------------------------------------------------------------------------------

# Number of time steps.
T = Vector{Int64}(1:nrow(df_weather_time))
# Set of electricity production facilities: dispatchable, storage and indispatchable facilities.
P = string.(df_supply_tech.technology)
# Set of already installed facilities for which the installed capacity will not be optimised, i.e., is fixed.
INSTP = string.(df_supply_tech[df_supply_tech.g_max .> 0, "technology"])
# Set of dispatchable electricity production facilities.
DISP = string.(df_supply_tech[df_supply_tech.dispatchable .== 1, "technology"])
# Set of non-dispatchable electricity production facilities.
# Contrary to the definition in the methodology, electricity storage facilities are contained in this set.
NDISP = string.(df_supply_tech[df_supply_tech.dispatchable .== 0, "technology"])
# Set of electricity storage facilities.
S = string.(df_supply_tech[df_supply_tech.inv_storage .> 0, "technology"])

# Get availability from the weather data DataFrame for wind and solar.
avail_types = Vector{String}(["pv", "wind"])
avail = DenseAxisArray{Float64}(undef, T, avail_types)
for i in T
    avail[i, "wind"] = coalesce(df_weather_time.windspeed[i],0)
    avail[i, "pv"] = coalesce(df_weather_time.solarenergy[i],0)
end
# Normalise the values for availability.
avail[:, "wind"] = avail[:, "wind"] ./ maximum(avail[:, "wind"])
avail[:, "pv"] = avail[:, "pv"] ./ maximum(avail[:, "pv"])

# Compute the marginal cost of electricity production of each dispatchable facility.
mc = map(eachrow(filter(x-> x.dispatchable == 1, df_supply_tech))) do row
    fc = row.fuel_cost
    ef = row.emission_factor
    eff = row.eff
    vc = row.variable_cost

    return fc/eff + co2_price*ef/eff + vc
end

mc = DenseAxisArray(mc, DISP)

# Save the efficiency factor of each electricity producing facility.
eff = DenseAxisArray(df_supply_tech.eff, P)

# Define two methods to compute the annualised investment cost per unit
# of capacity of electricity producing facilities.
annuity(i,lifetime) = i * ((1 + i)^lifetime) / (((1 + i)^lifetime) - 1)
calc_inv_cost(oc, i, lt) = oc * annuity(i, lt)

# Compute the annualised investment cost per unit of capacity of electricity producing facilities.
invest_cost = map(
    row-> row.g_max > 0 ? row.fixed_cost : calc_inv_cost(maximum([row.inv_power, row.inv_storage]),
                                                         interest_rate,
                                                         row.lifetime) + row.fixed_cost,
    eachrow(df_supply_tech)
)
invest_cost = DenseAxisArray(invest_cost, P)

# Maximum capacity of installed facilities.
inst_g_max = map(
    row -> row.g_max,
    eachrow(filter(x-> x.g_max > 0, df_supply_tech))
)
inst_g_max = DenseAxisArray(inst_g_max, INSTP)

# Surface covered per unit of capacity.
SW = DenseAxisArray(df_supply_tech.surface_covered, P)
# Maximum surface available for each electricity producing facility.
MSW = DenseAxisArray(df_supply_tech.max_surface, P)


# Set of water facilities.
W = string.(df_water_tech.technology)
# Set of water tanks.
WTa = string.(df_water_tech[df_water_tech.storage .== 1, "technology"])
# Set of water pumps.
WPu = string.(df_water_tech[df_water_tech.storage .== 0, "technology"])

# Cost of pumping one unit of water volume in units of electricity.
WEC = DenseAxisArray(df_water_tech.energy_cost, W)
# Water tank's capacities.
WTC = map(
    row -> row.max_capacity,
    eachrow(filter(x-> x.storage > 0, df_water_tech))
)
WTC = DenseAxisArray(WTC, WTa)

# Convert wind speed from km/h to m/s.
df_weather_time[!, :windspeed] = df_weather_time[!, :windspeed] ./ 3.6

# Delta: slope vapour pressure curve.
df_weather_time[!, :delta] =  map(
    row -> 4098 * 0.6108
           * exp(17.27*row.temp/(row.temp+273.3))
           / (row.temp+273.3)^2,
    eachrow(df_weather_time)
)

# R_{ns}: net solar radiation.
df_weather_time[!, :Rns] = map(
    row -> (1 - alpha) * row.solarradiation,
    eachrow(df_weather_time)
)

df_weather_time[!,:ea] = map(
    row -> 0.6108*exp(17.27*row.temp/(row.temp+237.3)) # e^0(T)
           * row.humidity/100,
    eachrow(df_weather_time)
)

# R_{nl}: net long wave radiation.
df_weather_time[!, :Rnl] = map(
    row -> 4.903 * 10^(-9)
            *(row.temp + 273.16)^4
            * (0.34 - 0.14 * sqrt(row.ea) * (1.35 * 0.67 - 0.35)),
    eachrow(df_weather_time)
)

# R_n: net radiation at the crop surface.
df_weather_time[!, :Rn] = map(
    row -> row.Rns - row.Rnl,
    eachrow(df_weather_time)
)

# G: soil heat flux density.
df_weather_time[!, :G] = map(
    row -> row.Rn * 0.1,
    eachrow(df_weather_time)
)

# Gamma: psychometric constant.
gamma = 0.665*10^(-3)* 101.3*((293 - 0.0065*altitude)/293)^(5.26)

# Hourly reference evapotranspiration in mm.
evapotranspiration = map(
    row-> (0.408*row.delta*(row.Rn - row.G)
           + gamma * 900/(row.temp+273) * row.windspeed * (0.6108*exp(17.27*row.temp/(row.temp+237.3))
           - row.ea))
    / (row.delta + gamma * (1 + 0.34*row.windspeed)),
    eachrow(df_weather_time)
)
minevotrans = minimum(evapotranspiration)
println("Min. evapotranspiration: $minevotrans.")

# Shift evapotranspiration figures to get all non-negative numbers.
if sum(evapotranspiration .< 0) > 0
    WE = evapotranspiration .+ abs.(minimum(evapotranspiration))
end
# Substract the precipitation in mm from the evapotranspiration.
WE = WE - df_weather_time.precip
# Normalise the irrigation efficiency.
WE = WE ./ maximum(WE)
# Compute the irrigation efficiency as the inverse of the evapotranspiration.
WE = 1 .- WE

# Add up the demand from all other demand sources hourly.
OD = sum(eachcol(df_other_demand[:, Cols(x -> startswith(x, "load"))]))

# -----------------------------------------------------------------------------------------------------------
# Model definitions                                                                                         |
# -----------------------------------------------------------------------------------------------------------

# Define a new JuMP model.
m = Model(GLPK.Optimizer)

# Define the model variables.
@variables m begin
    # Generation from dispatchable facilities (dispatchable generators AND storages).
    # Contrary to the definition in the methodology, G is not defined for non-dispatchable facilities.
    G[DISP, T] >= 0
    # Energy used for charging electricity storages.
    CH[S, T] >= 0
    # Electricity storage levels.
    L[S, T] >= 0
    # Installed maximum capacities.
    PW[P] >= 0
    # Water pumped by all water facilities.
    WP[W, T] >= 0
    # Water used for filling water tanks.
    WCH[WTa, T] >= 0
    # Water tank levels.
    WL[WTa, T] >= 0
    #
    # Variables to improve the probability of finding a feasible solution.
    # Curtailement.
    CU[T] >= 0
    # Lost load.
    LL[T] >= 0
end

# Expresion for the electricity production from non-dispatchable facilities.
@expression(m, feedin[p=NDISP, t=T],
    PW[p] * (startswith(p, "pv") ?  avail[t,"pv"] : avail[t,"wind"])
)

# Expression for the hourly electricity demand from water pumps and water tanks.
@expression(m, elec_demand_irrigation[t=T],
    sum(WP[w,t] * WEC[w] for w in W)
)

# Expression with the yearly electricity cost without investments.
@expression(m, yearly_elec_cost,
    8760/length(T) * (
        sum(mc[p] * G[p,t] for p in DISP , t in T)
        + sum(lost_load_cost * LL[t] for t in T)
    )
)

# Expression with the yearly investment cost.
@expression(m, yearly_investment_cost,
    sum(invest_cost[p] * PW[p] for p in P)
)

# Expression with the total water irrigated to the field over the optimisation period.
@expression(m, total_water_irrigated,
    sum(WE[t] * (sum(WP[w,t] for w in W)
        - sum(WCH[wta, t] for wta in WTa))
        for t in T)
)

# Define the objective function of the model.
@objective(m, Min,
    8760/length(T) * (
        sum(mc[p] * G[p,t] for p in DISP , t in T)
        + sum(lost_load_cost * LL[t] for t in T)
    )
    + sum(invest_cost[p] * PW[p] for p in P)
)

# Energy balance constraints.
@constraint(m,
    EnergyBalance[t=T],
    sum(G[p,t] for p in DISP)
    + sum(feedin[nd, t] for nd in NDISP)
    + LL[t]
    ==
    elec_demand_irrigation[t]
    + OD[t]
    + sum(CH[s,t] for s in S)
    + CU[t]
)

# Maximum generation constraints for dispatchable generators
@constraint(m, MaxGeneration[p=DISP,t=T], G[p,t] <= PW[p])
# Fixed maximum generation from already installed facilities.
@constraint(m, InstalledGeneration[inst=INSTP], PW[inst] == inst_g_max[inst])

# Maximum electricity storage constraints for charging.
@constraint(m, MaxCharge[s=S,t=T], CH[s,t] <= PW[s])
# Maximum electricity storage constraints for storage levels.
@constraint(m, MaxStorage[s=S,t=T], L[s,t] <= PW[s])

# Define a function that gives us the next hour of the set T
next_hour(t) = t == T[end] ? T[1] : T[findfirst(isequal(t), T) + 1]
# Electricity storage facility level difference constraints.
@constraint(m, StorageBalance[s=S, t=T],
    L[s,next_hour(t)]
    ==
    L[s,t]
    + eff[s] * CH[s,t]
    - G[s,t]/eff[s]
)

# Electricity storage initial level constraints.
@constraint(m, StartingStorageLevel[s=S],
    L[s,1] == 0
)

# Contraints for area covered by generation facilities.
@constraint(m, AreaCovered[p=P],
    PW[p] * SW[p] <= MSW[p]
)

# Constaints for minimum irrigation.
Days = Vector{Int64}(1:nrow(df_weather_time)/24)
get_day_hours(d) = Vector{Int64}(24*(d-1)+1 : 24*d)
@constraint(m, IrrigationBalance[d=Days],
    sum(WE[t] * (sum(WP[w,t] for w in W)
                 - sum(WCH[wta, t] for wta in WTa))
        for t in get_day_hours(d))
    >=
    DEW
)

# Maximum water storage constaints for charging.
@constraint(m, MaxWaterCharge[wta=WTa,t=T],
    WCH[wta,t] <= WTC[wta]
)
# Maximum water storage constraints for storage levels.
@constraint(m, MaxWaterStorage[wta=WTa,t=T],
    WL[wta,t] <= WTC[wta]    
)
# Water storage facility level difference constraints.
@constraint(m, WaterStorageBalance[wta=WTa,t=T],
    WL[wta,next_hour(t)]
    ==
    WE[t] * WL[wta,t] + WCH[wta,t] - WP[wta,t]
)

# Water storage initial level constraints.
@constraint(m, StartingWaterStorageLevel[wta=WTa],
    WL[wta,1] == 0
)

# -----------------------------------------------------------------------------------------------------------
# Optimise and print summary of the optmisation process and result.                                         |
# -----------------------------------------------------------------------------------------------------------

optimize!(m)
println(solution_summary(m))

# -----------------------------------------------------------------------------------------------------------
# Compute and print some statistical figures.                                                               |
# -----------------------------------------------------------------------------------------------------------

include(joinpath("Statistical_Figures.jl"))
Compute_Model_Statistics()

# -----------------------------------------------------------------------------------------------------------
# Plot results                                                                                              |
# -----------------------------------------------------------------------------------------------------------

include(joinpath("Plot.jl"))
fig = plot_result();