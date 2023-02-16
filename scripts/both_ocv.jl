#assume that both OCV's are already there.

x⁻₀ = 0.97
V⁺₀ = 4.2 + calcocv(anodeocv,x⁻₀,297.0)
x⁺₀ = OCV.get_x_from_voltage(cathodeocv,V⁺₀,297.0)
