using Revise

using Unitful, DifferentialEquations
#using Clapeyron
using Plots

using ModelingToolkit
using MethodOfLines
using DomainSets

R = ustrip(uconvert(u"(bar*m^3)/(mol*K)", 8.314u"J/(mol*K)"))
cross_sect_area = ustrip(uconvert(u"m^2", 10u"m^2"))
ρ_bulk = ustrip(uconvert(u"kg/m^3", 1000u"kg/m^3"))
Dp = ustrip(uconvert(u"m", 3u"mm"))
#NOTE: making the reactor larger in any dimension prevents t being set as a 1 dimensional array for some reason
epsilon = 0.3 #Bed void fraction
mu_gas = ustrip(uconvert(u"bar*s", 2.0e-5u"Pa*s"))

MW_CH3OH = ustrip(uconvert(u"kg/mol", 0.03204u"g/mol"))
MW_H2O = ustrip(uconvert(u"kg/mol", 0.018015u"g/mol"))
MW_CO = ustrip(uconvert(u"kg/mol", 0.02801u"g/mol"))
MW_H2 = ustrip(uconvert(u"kg/mol", 0.002016u"g/mol"))
MW_CO2 = ustrip(uconvert(u"kg/mol", 0.04401u"g/mol"))

kf_A_SMR = ustrip(uconvert(u"mol/(bar^2*kg*s)", 300000u"mol/(bar^2*kg*s)"))
Ea_f_SMR = ustrip(uconvert(u"(bar*m^3)/mol", 90u"kJ/mol"))
kr_A_SMR = ustrip(uconvert(u"mol/(bar^4*kg*s)", 0.008u"mol/(bar^4*kg*s)"))
Ea_r_SMR = ustrip(uconvert(u"(bar*m^3)/mol", 43u"kJ/mol"))
Keq_ref_SMR = 1000
ΔH_SMR = ustrip(uconvert(u"(bar*m^3)/mol", 49.4u"kJ/mol"))
ref_T_SMR = 298

kf_A_WGS = ustrip(uconvert(u"mol/(bar^2*kg*s)", 3000u"mol/(bar^2*kg*s)"))
Ea_f_WGS = ustrip(uconvert(u"(bar*m^3)/mol", 60u"kJ/mol"))
kr_A_WGS = ustrip(uconvert(u"mol/(bar^2*kg*s)", 35000u"mol/(bar^2*kg*s)"))
Ea_r_WGS = ustrip(uconvert(u"(bar*m^3)/mol", 100u"kJ/mol"))
Keq_ref_WGS= 36
ΔH_WGS = ustrip(uconvert(u"(bar*m^3)/mol", -41.1u"kJ/mol"))
ref_T_WGS = 298

#@parameters cross_sect_area ρ_bulk Dp epsilon mu_gas MW_CH3OH MW_H2O MW_CO MW_H2 MW_CO2 kf_A_SMR Ea_f_SMR kr_A_SMR Ea_r_SMR Keq_ref_SMR kf_A_WGS Ea_f_WGS kr_A_WGS Ea_r_WGS Keq_ref_WGS
@independent_variables t z 
@variables F_CH3OH(..) F_H2O(..) F_CO(..) F_H2(..) F_CO2(..) T(..) P(..)

Dt = Differential(t)
Dz = Differential(z)

F_total = F_CH3OH(t, z) + F_H2O(t, z) + F_CO(t, z) + F_H2(t, z) + F_CO2(t, z)

y_CH3OH = F_CH3OH(t, z) / F_total 
y_H2O = F_H2O(t, z) / F_total 
y_CO = F_CO(t, z) / F_total 
y_H2 = F_H2(t, z) / F_total 
y_CO2 = F_CO2(t, z) / F_total

p_CH3OH = y_CH3OH * P(t, z)
p_H2O = y_H2O * P(t, z)
p_CO = y_CO * P(t, z)
p_H2 = y_H2 * P(t, z)
p_CO2 = y_CO2 * P(t, z)

kf_SMR = kf_A_SMR * exp(-Ea_f_SMR / (R * T(t, z)))
Keq_SMR = Keq_ref_SMR * exp((-ΔH_SMR / R) * (1/T(t, z) - 1/ref_T_SMR))
kr_SMR = kf_SMR / 1

kf_WGS = kf_A_WGS * exp(-Ea_f_WGS / (R * T(t, z)))
Keq_WGS = Keq_ref_WGS * exp((-ΔH_WGS / R) * (1/T(t, z) - 1/ref_T_WGS))
kr_WGS = kf_WGS / 1

r_net_SMR = kf_SMR * (p_CH3OH * p_H2O) - kr_SMR * (p_CO * p_H2^3)
r_net_WGS = kf_WGS * (p_CO * p_H2O) - kr_WGS * (p_CO2 * p_H2)

source_CH3OH = ρ_bulk * cross_sect_area * (-1 * r_net_SMR)
source_H2O = ρ_bulk * cross_sect_area * (-1 * r_net_SMR - 1 * r_net_WGS)
source_CO = ρ_bulk * cross_sect_area * (1 * r_net_SMR - 1 * r_net_WGS)
source_H2 = ρ_bulk * cross_sect_area * (3 * r_net_SMR + 1 * r_net_WGS)
source_CO2 = ρ_bulk * cross_sect_area * (1 * r_net_WGS)

interstitial_velocity = (R * T(t, z) * F_total) / (P(t, z) * epsilon * cross_sect_area)

eq_F_CH3OH = Dt(F_CH3OH(t, z)) ~ interstitial_velocity * (-Dz(F_CH3OH(t, z)) + source_CH3OH)
eq_F_H2O = Dt(F_H2O(t, z)) ~ interstitial_velocity * (-Dz(F_H2O(t, z)) + source_H2O)
eq_F_CO = Dt(F_CO(t, z)) ~ interstitial_velocity * (-Dz(F_CO(t, z)) + source_CO)
eq_F_H2 = Dt(F_H2(t, z)) ~ interstitial_velocity * (-Dz(F_H2(t, z)) + source_H2)
eq_F_CO2 = Dt(F_CO2(t, z)) ~ interstitial_velocity * (-Dz(F_CO2(t, z)) + source_CO2)

MW_avg = y_CH3OH * MW_CH3OH + y_H2O * MW_H2O + y_CO * MW_CO + y_H2 * MW_H2 + y_CO2 * MW_CO2
mass_flow_total = F_CH3OH(t, z) * MW_CH3OH + F_H2O(t, z) * MW_H2O + F_CO(t, z) * MW_CO + F_H2(t, z) * MW_H2 + F_CO2(t, z) * MW_CO2
G = mass_flow_total / cross_sect_area
rho_gas = (P(t, z) * MW_avg) / (R * T(t, z))
ergun_term_1 = 150 * mu_gas * (1 - epsilon)^2 / (Dp^2 * epsilon^3)
ergun_term_2 = 1.75 * G * (1 - epsilon) / (Dp * epsilon^3)
ergun_expression = - (G / rho_gas) * (ergun_term_1 + ergun_term_2)

eq_P = Dz(P(t, z)) ~ ergun_expression

eq_T = Dt(T(t, z)) ~ -0.001

eqs = [eq_F_CH3OH, eq_F_H2O, eq_F_CO, eq_F_H2, eq_F_CO2, eq_P, eq_T]

t_start = 0
t_end = 200.0

z_min = 0
z_max = 1.0

domains = [t ∈ Interval(t_start, t_end),
           z ∈ Interval(z_min, z_max)]
bcs = [
    F_CH3OH(t_start, z) ~ 0.0001, #intial conditions throughout the z of the reactor
    F_H2O(t_start, z) ~ 0.0001,
    F_CO(t_start, z) ~ 0.0001,
    F_H2(t_start, z) ~ 0.0001,
    F_CO2(t_start, z) ~ 0.0001,
    P(t_start, z) ~ 1.0,
    T(t_start, z) ~ 573.0,
    F_CH3OH(t, z_min) ~ 100.0, #boundary conditions at the start of the reactor 
    F_H2O(t, z_min) ~ 100.0,
    F_CO(t, z_min) ~ 100.0,
    F_H2(t, z_min) ~ 100.0,
    F_CO2(t, z_min) ~ 100.0,
    P(t, z_min) ~ 1.0,
    T(t, z_min) ~ 573.0,
]
"""
Dz(F_CH3OH(t, z_min)) ~ 1,
    Dz(F_H2O(t, z_min)) ~ 0,
    Dz(F_CO(t, z_min)) ~ 0,
    Dz(F_H2(t, z_min)) ~ 0,
    Dz(F_CO2(t, z_min)) ~ 0,
    Dz(P(t, z_min)) ~ 0,
    Dz(T(t, z_min)) ~ 0,
Dz(F_CH3OH(t, z_max)) ~ 0,
    Dz(F_H2O(t, z_max)) ~ 0,
    Dz(F_CO(t, z_max)) ~ 0,
    Dz(F_H2(t, z_max)) ~ 0,
    Dz(F_CO2(t, z_max)) ~ 0,
    Dz(P(t, z_max)) ~ 0,
    Dz(T(t, z_max)) ~ 0,
"""

params = [R, cross_sect_area, ρ_bulk, Dp, epsilon, mu_gas, MW_CH3OH, 
          MW_H2O, MW_CO, MW_H2, MW_CO2, kf_A_SMR, Ea_f_SMR, kr_A_SMR, 
          Ea_r_SMR, Keq_ref_SMR, kf_A_WGS, Ea_f_WGS, kr_A_WGS, 
          Ea_r_WGS, Keq_ref_WGS] #Not necessary anymore 
@named pdesystem = PDESystem(eqs, bcs, domains, [t, z], [F_CH3OH(t, z), F_H2O(t, z), F_CO(t, z), F_H2(t, z), F_CO2(t, z), T(t, z), P(t, z)])

pdesystem.domain
pdesystem.bcs

discretization = MOLFiniteDifference([z => 10], t) 
@time prob = discretize(pdesystem, discretization)

@time sol = solve(prob, Rosenbrock23())

#Plot results
z_grid = sol[z]
t_grid = sol[t]
F_CH3OH_sol = sol[F_CH3OH(t, z)]
P_sol = sol[P(t, z)]
T_sol = sol[T(t, z)]
# etc. for other variables

# Plot pressure drop
plt1 = plot(z_grid, sol[P(t, z)][1, :], xlabel="Reactor Length (m)", ylabel="Pressure (bar)", label="Pressure", lw=2)

plot(t_grid, sol[F_CH3OH(t, z)][:, 1], xlabel="time (s)", ylabel="Molar Flow (mol/s)", label="CH3OH", lw=2)
sol[F_CH3OH(t, z)][:, 2]
plot(t_grid, sol[F_CH3OH(t, z)][:, 1], xlabel="time (s)", ylabel="Molar Flow (mol/s)", label="CH3OH", lw=2)
sol[F_H2(t, z)][:, 2]

# Plot molar flows
time_for_plot = 18
current_time = t_grid[time_for_plot]
z_slice = z_grid[4:end]
plot(z_slice, sol[F_CH3OH(t, z)][time_for_plot, 4:end], xlabel="Reactor Length (m)", ylabel="Molar Flow (mol/s)", label="CH3OH", lw=2)
sol[F_CH3OH(t, z)][time_for_plot, 4:end]
plot!(z_slice, sol[F_H2O(t, z)][time_for_plot, 4:end], label="H2O")
sol[F_H2O(t, z)][time_for_plot, 4:end]
plot!(z_slice, sol[F_H2(t, z)][time_for_plot, 4:end], label="H2")
sol[F_H2(t, z)][time_for_plot, 4:end]
plot!(z_slice, sol[F_CO(t, z)][time_for_plot, 4:end], label="CO")
sol[F_CO(t, z)][time_for_plot, 4:end]
plot!(z_slice, sol[F_CO2(t, z)][time_for_plot, 4:end], label="CO2")
sol[F_CO2(t, z)][time_for_plot, 4:end]