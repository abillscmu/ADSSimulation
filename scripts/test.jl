using JLD2
using ADSSimulation
using PythonCall
using DataFrames
using CSV

pybamm = pyimport("pybamm")

@load "/Users/abills/Datasets/ADS_CMU_Ageing/CMU/room_temperature.jld2"

key = Dict(
    "First Cycle" => 1,
    "Serial Number" => 1,
    "File Index" => "01"
)

df = room_temperature[key]

cla()
plot(df[!,"time/s"],df[!,"I/mA"])