using DataFrames, CSV, OCV, PythonPlot, JuMP
#pygui(true)


sic_df = CSV.read("/Users/abills/Downloads/SiC.csv", DataFrame, header=0)

c_s_n_max = 34684.
c_s_n_min = 0.0

x_0_100 = 0.852 * (sic_df.Column1) .- minimum(sic_df.Column1)

neg_x_spline = LinearInterpolation(x_0_100, sic_df.Column1)

anodeocv = OCV.SplineOCV(x_0_100, sic_df.Column2, AkimaInterpolation, c_s_n_max, c_s_n_min)


figure(1)
clf()
scatter(x_0_100, sic_df.Column2)
V = anodeocv.(x_0_100)
plot(x_0_100,V)
grid()
legend(["Experiment(SiC)","Akima Spline"])
xlabel("Stoichiometry")
ylabel("Negative Electrode Potential [V]")


