import pybamm

def negative_ocp(sto):
        u_eq = (
            1.9793 * pybamm.exp(-39.3631 * sto)
            + 0.2482
            - 0.0909 * pybamm.tanh(29.8538 * (sto - 0.1234))
            - 0.04478 * pybamm.tanh(14.9159 * (sto - 0.2769))
            - 0.0205 * pybamm.tanh(30.4444 * (sto - 0.6103))
        )

        return u_eq

def negative_exchange_current_density(c_e, c_s_surf, c_s_max, T):
    m_ref = 1e-11 * pybamm.constants.F
    E_r = 3600
    arrhenius = pybamm.exp(E_r * (1 / 298.15 - 1 / T))
    return (
        m_ref
        * arrhenius
        * c_e**0.5
        * c_s_surf**0.5
        * (c_s_max - c_s_surf) ** 0.5
    )

def positive_ocp(sto):
    u_eq = (
        -0.8090 * sto
        + 4.4875
        - 0.0428 * pybamm.tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * pybamm.tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * pybamm.tanh(15.9308 * (sto - 0.3120))
    )
    return u_eq

def positive_exchange_current_density(c_e, c_s_surf, c_s_max, T):
    m_ref = 3e-11 * pybamm.constants.F
    E_r = 3600
    arrhenius = pybamm.exp(E_r * (1 / 298.15 - 1 / T))
    return (
        m_ref
        * arrhenius
        * c_e**0.5
        * c_s_surf**0.5
        * (c_s_max - c_s_surf) ** 0.5
    )

def electrolyte_diffusivity_Valoen2005(c_e, T):
    # mol/m3 to molar
    c_e = c_e / 1000
    T_g = 229 + 5 * c_e
    D_0 = -4.43 - 54 / (T - T_g)
    D_1 = -0.22
    # cm2/s to m2/s
    # note, in the Valoen paper, ln means log10, so its inverse is 10^x
    return (10 ** (D_0 + D_1 * c_e)) * 1e-4

def electrolyte_conductivity_Valoen2005(c_e, T):
    # mol/m3 to molar
    c_e = c_e / 1000
    # mS/cm to S/m
    return (1e-3 / 1e-2) * (
        c_e
        * (
            (-10.5 + 0.0740 * T - 6.96e-5 * T**2)
            + c_e * (0.668 - 0.0178 * T + 2.80e-5 * T**2)
            + c_e**2 * (0.494 - 8.86e-4 * T)
        )
        ** 2
    )
    
x_100 = 0.852
c_n_max = 34684
y_100 = 0.222
c_p_max = 50060
parameters = {
    # negative electrode
    "Maximum concentration in negative electrode [mol.m-3]": c_n_max,
    "Negative electrode thickness [m]": 86.7e-6,
    "Negative electrode active material volume fraction": 0.694,
    "Negative electrode porosity": 0.216,
    "Negative electrode OCP [V]": negative_ocp,
    "Negative electrode exchange-current density [A.m-2]": negative_exchange_current_density,
    "Negative electrode OCP entropic change [V.K-1]": 0,
    "Initial concentration in negative electrode [mol.m-3]": c_n_max * x_100,
    "Negative particle radius [m]": 6.1e-6,
    "Negative electrode Bruggeman coefficient (electrolyte)": 1.5,
    "Negative electrode Bruggeman coefficient (electrode)": 1.5,
    "Negative electrode conductivity [S.m-1]": 100,
    "Negative electrode diffusivity [m2.s-1]": 5e-14,
    # Separator
    "Separator thickness [m]": 12e-6,
    "Separator porosity": 0.45,
    "Separator Bruggeman coefficient (electrolyte)": 1.5,
    # Positive electrode
    "Maximum concentration in positive electrode [mol.m-3]": c_p_max,
    "Positive electrode thickness [m]": 66.2e-6,
    "Positive electrode active material volume fraction": 0.745,
    "Positive electrode porosity": 0.171,
    "Positive electrode OCP [V]": positive_ocp,
    "Positive electrode exchange-current density [A.m-2]": positive_exchange_current_density,
    "Positive electrode OCP entropic change [V.K-1]": 0,
    "Initial concentration in positive electrode [mol.m-3]": c_p_max * y_100,
    "Positive particle radius [m]": 3.8e-6,
    "Positive electrode Bruggeman coefficient (electrolyte)": 1.85,
    "Positive electrode Bruggeman coefficient (electrode)": 1.85,
    "Positive electrode conductivity [S.m-1]": 0.17,
    "Positive electrode diffusivity [m2.s-1]": 5e-13,
    # electrolyte
    "Typical electrolyte concentration [mol.m-3]": 1000.0,
    "Initial concentration in electrolyte [mol.m-3]": 1000.0,
    "Cation transference number": 0.38,
    "Electrolyte diffusivity [m2.s-1]": electrolyte_diffusivity_Valoen2005,
    "Electrolyte conductivity [S.m-1]": electrolyte_conductivity_Valoen2005,
    "1 + dlnf/dlnc": 1.0,
    # cell
    "Electrode height [m]": 5.8e-2,
    "Electrode width [m]": 61.5e-2 * 2,
    "Lower voltage cut-off [V]": 2.5,
    "Upper voltage cut-off [V]": 4.2,
    "Typical current [A]": 3.35,
    "Current function [A]": 3.35,
    "Nominal cell capacity [A.h]": 3.35,
    "Initial concentration in negative electrode [mol.m-3]": 30163.186289632762,
    "Maximum stoichiometry in negative electrode": 0.8696570836591155,
    "Minimum stoichiometry in negative electrode": 0.030109799223134615,
    "Initial concentration in positive electrode [mol.m-3]": 13208.233442104221,
    "Minimum stoichiometry in positive electrode": 0.2638480511806676,
    "Maximum stoichiometry in positive electrode": 0.9735039208548855
}

parameter_values = pybamm.ParameterValues("Chen2020")
parameter_values.update(parameters, check_already_exists=False)