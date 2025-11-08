using Unitful, DifferentialEquations
using Plots

using ModelingToolkit


@variables W

@connector function FlowPort(; species::Vector{Symbol}, name) 
    @parameters begin
        species = [], [description="the order of the components must match the stoichiometric coefficients inputted in ReactionKinetics, format would be [CH3OH, ]"]
    end
    @variables begin
        T(W)
        P(W)
        (F(W))[species]
    end

    for s in species
        setmetadata(vars.F[s], :connect, Flow)
    end
    return ODESystem(Equation[], W; name=name, states=[vars.T, vars.P, vars.F...])
end

@mtkmodel ReactionKinetics begin
    @parameters begin
        Ea_f
        Ea_r
        kf_A
        kf_Ea
        kr_A
        kr_Ea
        reactants_stoichs_vec, [description="vector matching "]
        products_stoichs_vec, [description="Overall heat transfer coeff"]
        R_gas = 8.314
    end
    
    @variables begin
        T(W)
        P_total(W)
        rate(W)
        reactants_partial_pressures(W)[1:length(reactants_stoichs_vec)]
        products_partial_pressures(W)[1:length(products_stoichs_vec)]
    end
    
    @equations  begin
        kf = kf_A * exp(-kf_Ea / (R_GAS * T))
        kr = kr_A * exp(-kr_Ea / (R_GAS * T))
        forward_term[1:length(reactants_stoichs_vec)] *= [reactants_partial_pressures[reactnant]^reactants_stoichs_vec[reactant] for reactant in eachindex(reactants_stoichs_vec)]
        forward_term[1:length(products_stoichs_vec)] *= [products_partial_pressures[product]^products_stoichs_vec[product] for product in eachindex(product_stoichs_vec)]
        rate ~ kf .* forward_term - kr .* reverse_term
    end
end

@mtkmodel MethanolDecompositionKinetics begin
    @parameters begin
        kf_A = 1.39e5,  [description="Forward pre-exponential factor"]
        kf_Ea = 90000.0, [description="Forward activation energy J/mol"]
        kr_A = 1.0      # Set to a non-zero value to avoid issues if p_CO or p_H2 is zero initially
        kr_Ea = 43000.0, [description="Reverse activation energy J/mol"]
        R_GAS = 8.314,   [description="Universal gas constant"]
    end
    
    # These are the inputs this component needs from the reactor environment
    @variables begin
        T(W),        [description="Local Temperature"]
        p_CH3OH(W),  [description="Partial pressure of Methanol"]
        p_CO(W),     [description="Partial pressure of CO"]
        p_H2(W),     [description="Partial pressure of H2"]
        
        # This is the output of the component
        rate(W),     [description="Net reaction rate mol/kg_cat/s"]
    end
    
    @equations  begin
        # Arrhenius equations
        kf = kf_A * exp(-kf_Ea / (R_GAS * T))
        kr = kr_A * exp(-kr_Ea / (R_GAS * T))
        
        # Rate Law for CH3OH <=> CO + 2*H2
        rate ~ kf * p_CH3OH - kr * p_CO * p_H2^2
    end
end

@mtkmodel WaterGasShiftKinetics begin
    @parameters begin
        kf_A = 1.02e-1
        kf_Ea = 60000.0
        kr_A = 2.55e4
        kr_Ea = 100000.0
        R_GAS = 8.314
    end
    
    @variables begin
        T(W) 
        p_CO(W)
        p_H2O(W)
        p_CO2(W)
        p_H2(W)
        rate(W)
    end
    
    @equations begin
        # Arrhenius equations
        kf = kf_A * exp(-kf_Ea / (R_GAS * T))
        kr = kr_A * exp(-kr_Ea / (R_GAS * T))
        
        # Rate Law for CO + H2O <=> CO2 + H2
        rate ~ kf * p_CO * p_H2O - kr * p_CO2 * p_H2
    end
end

function symbolic_ergun_pressure_drop(P, T, F_total, params)
    # This is a placeholder for the full symbolic Ergun equation.
    # For simplicity, let's just make it proportional to flow and inversely to P.
    G = params.G_const * F_total # Superficial velocity approximation
    beta0 = G * (150 * (1 - params.eps) * params.mu / params.Dp^2 + 1.75 * G / params.Dp)
    return -beta0 / (params.rho_cat * P) # Simplified form
end


@mtkmodel PackedBedReactorSegment begin
    @parameters begin
        # Heat transfer params
        U = 0.0,       [description="Overall heat transfer coeff W/m^2/K"]
        a = 0.0,       [description="Reactor surface area per kg catalyst m^2/kg"]
        T_amb = 298.15, [description="Ambient temperature K"]
        Q_in = 0.0,      [description="External heat input per kg W/kg"]
        
        # Pressure drop params
        G_const = 0.002, [description="Factor to relate total flow to mass velocity"]
        eps = 0.3,     [description="Bed void fraction"]
        mu = 2.0e-5,   [description="Viscosity Pa*s"]
        Dp = 0.005,    [description="Particle diameter m"]
        rho_cat = 1000.0,[description="Catalyst density kg/m^3"]

        # Heats of reaction J/mol
        H_rxn_MD = 90700.0
        H_rxn_WGS = -41100.0

        species = ["CH3OH", "H2O", "CO", "CO2", "H2"]
    end

    @components begin
        # Interface ports, now correctly initialized with the species list
        inlet = FlowPort(; species)
        outlet = FlowPort(; species)

        # Instantiate the reaction kinetics for each reaction
        MD_rxn = MethanolDecompositionKinetics()
        WGS_rxn = WaterGasShiftKinetics()
    end
    
    # These are the internal states of the reactor, varying with W
    # The F(W)[species] syntax creates F_CH3OH, F_H2O, etc. automatically
    @variables begin
        T(W)
        P(W)
        F(W)[1:length(species)]
    end

    @equations begin
        # --- Intermediate Variables ---
        F_total = sum(F)
        
        # --- Connect Reactor State to Kinetics Inputs ---
        # Connect temperature to both reactions
        MD_rxn.T ~ T
        WGS_rxn.T ~ T
        
        # Connect partial pressures to each reaction
        # p_i = y_i * P = (F_i / F_total) * P
        MD_rxn.p_CH3OH ~ (F[:CH3OH] / F_total) * P
        MD_rxn.p_CO     ~ (F[:CO]     / F_total) * P
        MD_rxn.p_H2     ~ (F[:H2]     / F_total) * P
        
        WGS_rxn.p_CO   ~ (F[:CO]   / F_total) * P
        WGS_rxn.p_H2O  ~ (F[:H2O]  / F_total) * P
        WGS_rxn.p_CO2  ~ (F[:CO2]  / F_total) * P
        WGS_rxn.p_H2   ~ (F[:H2]   / F_total) * P

        # --- Calculate Net Rate of Formation for Each Species ---
        # r_i = sum(nu_ij * R_j) over all reactions j
        r_CH3OH = -1 * MD_rxn.rate
        r_H2O   = -1 * WGS_rxn.rate
        r_CO    = +1 * MD_rxn.rate - 1 * WGS_rxn.rate
        r_H2    = +2 * MD_rxn.rate + 1 * WGS_rxn.rate
        r_CO2   = +1 * WGS_rxn.rate
        
        # --- Core Differential Equations ---
        # Mass Balances: dF_i / dW = r_i
        D(F[:CH3OH]) ~ r_CH3OH
        D(F[:H2O])   ~ r_H2O
        D(F[:CO])    ~ r_CO
        D(F[:H2])    ~ r_H2
        D(F[:CO2])   ~ r_CO2

        # Energy Balance: dT / dW = ...
        # Simplified Cp (could be replaced with a symbolic function of T)
        # Cp units: J/mol/K
        Cp_vals = [48.0, 34.0, 29.0, 29.0, 37.0] # CH3OH, H2O, CO, H2, CO2
        Cp_total_flow = sum(F[s] * Cp for (s, Cp) in zip(species, Cp_vals))
        
        heat_from_reactions = (-H_rxn_MD * MD_rxn.rate) + (-H_rxn_WGS * WGS_rxn.rate)
        heat_exchange = U * a * (T_amb - T)
        D(T) ~ (heat_from_reactions + heat_exchange + Q_in) / Cp_total_flow
        
        # Momentum Balance: dP / dW = ...
        D(P) ~ symbolic_ergun_pressure_drop(P, T, F_total; Dp, eps, mu, rho_cat, G_const)

        # --- Outlet Port Connections ---
        # The outlet state is the same as the internal reactor state
        outlet.T ~ T
        outlet.P ~ P
        # Loop to connect all species flows
        [outlet.F[s] ~ F[s] for s in species]...
    end
end

# Define simulation parameters
W_span = (0.0, 100.0) # Integrate from 0 to 100 kg of catalyst

# Inlet conditions
inlet_T = 573.15  # K
inlet_P = 1.0e5   # Pa
inlet_flows = Dict(
    :CH3OH => 0.065,
    :H2O   => 0.065 * 1.3,
    :CO    => 1e-6,
    :H2    => 1e-6,
    :CO2   => 1e-6
)

# Instantiate the model
@mtkbuild reactor = PackedBedReactorSegment()

# Create the initial conditions vector u0
# The order must match the states of the system.
# MTK provides helpers to ensure this is correct.
species_list = reactor.species
u0 = [
    reactor.T => inlet_T,
    reactor.P => inlet_P,
    # Map each species to its inlet flow
    [reactor.F[s] => inlet_flows[s] for s in species_list]...
]

# Create the parameter map (if you want to override defaults)
p_map = [
    reactor.U => 10.0, # W/m^2/K
    reactor.a => 0.5,  # m^2/kg
    reactor.T_amb => 300.0 # K
]

# Create and solve the ODE problem
prob = ODEProblem(reactor, u0, W_span, p_map)
sol = solve(prob, Rodas5P())

# --- Plotting Results ---
# The solution object can be indexed by the symbolic variables!
plot(sol, idxs=[reactor.T], title="Temperature Profile", xlabel="Catalyst Weight (kg)", ylabel="Temperature (K)")
plot(sol, idxs=[reactor.P], title="Pressure Profile", xlabel="Catalyst Weight (kg)", ylabel="Pressure (Pa)")

# Plot molar flows
plot(sol, idxs=[reactor.F[s] for s in species_list], 
     title="Molar Flow Profiles", xlabel="Catalyst Weight (kg)", ylabel="Molar Flow (mol/s)",
     labels=permutedims(string.(species_list)))