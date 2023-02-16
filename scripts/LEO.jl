using PythonCall
using PythonPlot
using DataFrames
using CSV
using JLD2
using ADSSimulation
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
df = df[378:end,:]
df[!,"time/s"] = df[!,"time/s"] .- df[1, "time/s"]


sys.path.append(joinpath(@__DIR__,"../pysrc"))
pybamm_params = pyimport("pybamm_params")

parameter_values = pybamm_params.parameter_values
parameter_values["Current function [A]"] = 0.679


parameter_values["Maximum concentration in negative electrode [mol.m-3]"] = 31684
parameter_values["Maximum concentration in positive electrode [mol.m-3]"] = 48560



params = ADSSimulation.RoomTemperature["CYCLE_PARAMETERS"][RoomTemperature["CELLS"]["SN-01"]]

experiment = pybamm.Experiment(repeat([
    "Discharge at $(params["DISCHARGE_RATE"])mA $(params["EOD_CRITERIA"])",
    "Charge at $(params["CHARGE_RATE"])mA until $(params["VEOC"])V",
    "Hold at $(params["VEOC"])V until 170mA",
], 13))


spme = pybamm.lithium_ion.DFN()
spme_sim = pybamm.Simulation(spme, parameter_values=parameter_values, experiment=experiment)
sol = spme_sim.solve(initial_soc = 0.93)


figure(1)
clf()
plot(sol["Time [s]"].entries,sol["Terminal voltage [V]"].entries)
plot(df[!,"time/s"], df[!, "Ecell/V"])
legend(["Model","Experiment"])

xlim([-200,50000])