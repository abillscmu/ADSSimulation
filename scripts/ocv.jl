using PythonCall
using PythonPlot
using DataFrames
using CSV
using JLD2
using ADSSimulation


RAW_DATA_PATH = ADSSimulation.RAW_DATA_PATH

data = load(RAW_DATA_PATH * "/room_temperature.jld2")
room_temperature = data["room_temperature"]

serial_number = 1
first_cycle = 1
file_index = "04"

key = Dict(
    "First Cycle" => 1,
    "Serial Number" => 1,
    "File Index" => "01"
)

df = room_temperature[key]

df = filter(:Ns => ns -> ns >=8, df)
df[!,"time/s"] = df[!, "time/s"] .- df[1, "time/s"]

pygui(true)

pybamm = pyimport("pybamm")
sys = pyimport("sys")
sys.path.append(joinpath(@__DIR__,"../pysrc"))
pybamm_params = pyimport("pybamm_params")

parameter_values = pybamm_params.parameter_values
parameter_values["Current function [A]"] = 0.679


parameter_values["Maximum concentration in negative electrode [mol.m-3]"] = 34684
parameter_values["Maximum concentration in positive electrode [mol.m-3]"] = 50060



spm = pybamm.lithium_ion.SPM()
spme = pybamm.lithium_ion.SPMe()
dfn = pybamm.lithium_ion.DFN()

spm_sim = pybamm.Simulation(spm, parameter_values=parameter_values)
spme_sim = pybamm.Simulation(spme, parameter_values=parameter_values)
dfn_sim = pybamm.Simulation(spm, parameter_values=parameter_values)

spm_sim.build(initial_soc=1)
spme_sim.build(initial_soc=1)
dfn_sim.build(initial_soc=1)

spm_sol = spm_sim.solve(collect(0:20009))
spme_sol = spme_sim.solve(collect(0:20009))
dfn_sol = dfn_sim.solve(collect(0:20009))

figure(1)
clf()
plot(spm_sol["Time [s]"].entries,spm_sol["Terminal voltage [V]"].entries)
plot(spme_sol["Time [s]"].entries,spme_sol["Terminal voltage [V]"].entries)
plot(dfn_sol["Time [s]"].entries,dfn_sol["Terminal voltage [V]"].entries)
plot(df[!,"time/s"], df[!, "Ecell/V"],"k--")
xlabel("Time [s]")
ylabel("Terminal voltage [V]")
legend(["SPM", "SPMe","DFN","Experiment: SN-1, Cy-1"])
grid()
