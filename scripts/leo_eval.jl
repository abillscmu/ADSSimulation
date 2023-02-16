using PythonCall
using PythonPlot
using DataFrames
using CSV
using JLD2
using ADSSimulation
using CellFitElectrolyte
using OCV
using ComponentArrays
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using Statistics


cache = CellFitElectrolyte.initialize_cache(Float64)
@load "data/cathodeocv.jld2"
@load "data/anodeocv.jld2"

RoomTemperature = ADSSimulation.RoomTemperature
sys = pyimport("sys")
pybamm = pyimport("pybamm")
RAW_DATA_PATH = ADSSimulation.RAW_DATA_PATH
data = load(RAW_DATA_PATH * "/room_temperature.jld2")
room_temperature = data["room_temperature"]

key = Dict(
    "First Cycle" => 1,
    "Serial Number" => 1,
    "File Index" => "02"
)

df = room_temperature[key]
df = df[378:500,:]
df[!,"time/s"] = df[!,"time/s"] .- df[1, "time/s"]





initialcond = Dict("Starting Voltage[V]"=>4.1,"Ambient Temperature[K]" => 293.15)
current = -df[!,"I/mA"]./1000

current_interpolant = LinearInterpolation(current,df[!,"time/s"])
voltage_interpolant = LinearInterpolation(df[!,"Ecell/V"],df[!,"time/s"])

interpolated_time = collect(range(df[1,"time/s"],stop=df[end,"time/s"],step=1.0))


interpolated_current = current_interpolant.(interpolated_time)
interpolated_voltage = voltage_interpolant.(interpolated_time)


#set up cycle arrays
cycle_array = CellFitElectrolyte.current_profile(interpolated_current,interpolated_time)
num_steps = Int(cycle_array[1])
times = cycle_array[2:num_steps+1]
types = cycle_array[num_steps+2:num_steps+1+num_steps]
values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]



cellgeometry = ADSSimulation.cell_geometry()





function evaluator(p::ComponentVector{T}) where {T}
    # Handle Initial Conditions
    u::Array{T,1}  = Array{T,1}(undef,7)
    CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    Temp = 298
    @pack! p = Temp
    du = similar(u)
    input_type::T = types[1]
    input_value::T = values[1]
    @pack! p = input_type,input_value

    #Create Function and Initialize Integrator
    func = ODEFunction((du, u, p, t)->CellFitElectrolyte.equations_electrolyte_allocating(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv))
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false, tstops = times)

    #we're really only interested in temperature and voltage, so we'll just save those
    endV::Array{T,1} = Array{T,1}(undef,length(interpolated_voltage)-1)
    endt::Array{T,1} = Array{T,1}(undef,length(interpolated_voltage)-1)

    for step::Int in 1:num_steps-1
        Temp = 298.
        @pack! p = Temp
        input_type = types[step]
        input_value = values[step]
        end_time::T = times[step+1]
        @pack! p = input_type,input_value
        while integrator.t < times[step+1]
            step!(integrator)
        end
        if any(integrator.u .< 0)
            @warn "something is less than 0"
            return integrator
        end
        cₛˢ⁺ = integrator.u[7]
        cₛˢ⁻ = integrator.u[1]
        x⁺ = (cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min)
        x⁻ = (cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min)
        if (x⁺ >= 1)
            return integrator
        elseif (x⁻ >= 1)
            return integrator
        end
        endV[step] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
        endt[step] = integrator.t
    end
    return endV
end


function fit_cfe_ads(arr)
    frac_sol_am_neg = arr[1]
    frac_sol_am_pos = arr[2]
    εₑ⁺ = arr[3]
    εₑ⁻ = arr[4]
    x⁻₀ = arr[5]
    ω = arr[6]

    εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
    εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
    εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
    εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

    p = ComponentVector(θₛ⁻ = 1e-8, θₑ = 5e-7, θₛ⁺ = 6.547741580032837e-8, R⁺ = 3.8e-6, R⁻ = 6.1e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = x⁻₀, εₑˢ = 0.8, cₑ₀ = 1000.0, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)

    try
        V = evaluator(p)
        rmse = sqrt(mean((V .- interpolated_voltage[1:end-1]).^2))
        return V
    catch
        return 1e12
    end
end


frac_sol_am_neg = .885
frac_sol_am_pos = .9
εₑ⁺ = 0.171
εₑ⁻ = 0.216
x⁻₀ = 0.856
ω = 0.015

parent = [frac_sol_am_neg, frac_sol_am_pos, εₑ⁺, εₑ⁻, x⁻₀, ω]
ub = [1,1,0.5,0.5, 0.99, 0.05]
lb = [0.01, 0.01, 0.01, 0.01, 0.8, 0.0]

#CellFitElectrolyte.anneal(fit_cfe_ads, parent, ub, lb)



#println(rmse)

arr = [0.20678103076466636, 0.9015232723068044, 0.09767537701626247, 0.15770519030456417, 0.8303459498804079, 0.014800200150565782]

V = fit_cfe_ads(arr)


figure(1)
clf()
plot(interpolated_time[1:end-1], V)
plot(df[!,"time/s"],df[!,"Ecell/V"])
legend(["Model","Experiment"])

