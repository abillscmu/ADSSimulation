using DataFrames, CSV, OCV, PythonPlot, JuMP

c_s_p_max = 50060.
c_s_p_min = 0.0


nmc_df = CSV.read("/Users/abills/Downloads/NMC811_new.csv", DataFrame, header=0)

x_0_100 = (0.222 - 0.942).*nmc_df.Column1 .+ 0.942
#x_0_100 = 1 .- x_0_100

sorted_og_x = sort(nmc_df.Column1)

sorted_indices = sortperm(x_0_100)
sorted_x = x_0_100[sorted_indices]
sorted_V = nmc_df.Column2[sorted_indices]


cathodeocv = OCV.SplineOCV(sorted_x, sorted_V, AkimaInterpolation, c_s_p_max, c_s_p_min)

pos_x_spline = LinearInterpolation(sorted_x, sorted_og_x)

figure(1)
clf()
scatter(x_0_100, nmc_df.Column2)
V = cathodeocv.(sort(x_0_100))
plot(sort(x_0_100),V)
grid()
legend(["Experiment(NMC811)","Akima Spline"])
xlabel("Stoichiometry")
ylabel("Positive Electrode Potential [V]")

