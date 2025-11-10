using ModelingToolkit
using DifferentialEquations
#using Plots
#using GLMakie


import ModelingToolkit: t_nounits as t, D_nounits as D
#Start of Component Definitions
@connector ReactorPort begin
    @structural_parameters begin
        n_species
        n_reactants
        n_products
    end

    @parameters begin
        F_guess[1:n_species]
        T_guess = 273.13, [description = "temperature in K"]
        p_guess = 101325, [description = "pressure in Pa"]
    end

    @variables begin
        F(t)[1:n_species]#, [connect = Flow, guess = F_guess]
        T(t), [guess = T_guess]
        p(t), [guess = p_guess]
    end
end

@mtkmodel ReactionKinetics begin
    @structural_parameters begin
        n_species
        n_reactants
        n_products
        #reactant_indices[1:n_reactants]
        #product_indices[1:n_products]
    end

    @parameters begin
        reactants_stoichs[1:n_reactants]
        products_stoichs[1:n_products]
        kf_A, [description="Forward Arrhenius pre-exponential factor"]
        kf_Ea, [description="Forward reaction activation energy"]
        kr_A, [description="Reverse Arrhenius pre-exponential factor"]
        kr_Ea, [description="Reverse reaction activation energy"]
        R_gas = 8.314
    end
    
    @variables begin
        F(t)[1:n_species]
        T(t)
        p(t)
        rate(t)
    end
    
    @equations  begin
        F_total = sum(F)
        
        partial_pressures = p .* F ./ F_total
        
        kf = kf_A * exp(-kf_Ea / (R_gas * T))
        kr = kr_A * exp(-kr_Ea / (R_gas * T))

        forward_term = prod(partial_pressures[1]^reactants_stoichs[i] for i in 1:n_reactants)
        reverse_term = prod(partial_pressures[1]^products_stoichs[i] for i in 1:n_products)

        rate ~ kf * forward_term - kr * reverse_term
    end
end

@mtkmodel ErgunPressureDrop begin
    @parameters begin
        G = 0.002, [description="Factor to relate total flow to mass velocity"]
        ϵ = 0.3,     [description="Bed void fraction"]
        μ = 2.0e-5,   [description="Viscosity Pa*s"]
        Dp = 0.005,    [description="Particle diameter m"]
        ρ = 1000.0,[description="Catalyst density kg/m^3"]
    end
    
    @variables begin
        p(t)
    end

    @equations begin
        -D(p) ~ (150 * μ * G*(1 - ϵ)^2) / (Dp^2 * ϵ^3) + (1.75 * G^2 * (1 - ϵ)) / (Dp * ϵ^3)
    end
end

#Note sure if a reactor segment requires a 1D Element to connect them 
#Perhaps this would be a good use case for trying out the finite volume method in mtkcomponents
@mtkmodel ReactorSegment begin
    @components begin
        inlet = ReactorPort(n_species=n_species, n_reactants=n_reactants, n_products=n_products)
        outlet = ReactorPort(n_species=n_species, n_reactants=n_reactants, n_products=n_products)
        kinetics = ReactionKinetics(n_species=n_species, n_reactants=n_reactants, n_products=n_products)
        pressure_drop = ErgunPressureDrop()
    end

    @structural_parameters begin
        n_species
        n_reactants
        n_products
    end

    @parameters begin
        catalyst_weight = 1.0, [description = "weight of catalyst in kg present in reactor segment"] 
        heat_of_reaction = -92.3
        heat_capacity_vec[1:n_species]
        stoich_coeffs[1:n_species]
    end

    @variables begin
        F(t)[1:n_species], [description = "test"]
        T(t) = 293.13
        p(t) = 101325
    end

    @equations begin
        kinetics.F ~ F
        kinetics.T ~ T
        kinetics.p ~ p
        pressure_drop.p ~ p

        D.(F) ~ kinetics.rate .* stoich_coeffs

        D(T) ~ (kinetics.rate * -heat_of_reaction) / sum(F .* heat_capacity_vec)

        inlet.F ~ F
        inlet.T ~ T
        inlet.p ~ p
        
        outlet.F ~ F
        outlet.T ~ T
        outlet.p ~ p
    end
end


@mtkmodel SpeciesSource begin
    @components begin
        port = ReactorPort(n_species=n_species, n_reactants=n_reactants, n_products=n_products)
    end

    @structural_parameters begin
        n_species
        n_reactants
        n_products
    end

    @parameters begin
        F_fixed[1:n_species], [description = "molar flows in mol/s for each respective species"]
        T_fixed = 273.13, [description = "temperature in K"]
        p_fixed = 101325, [description = "pressure in Pa"]
    end

    @equations begin
        port.F ~ F_fixed
        port.T ~ T_fixed
        port.p ~ p_fixed
    end
end

stoich_coeffs_val = [-1, 0, 1, 2, 0] # Negative for reactants, positive for products
reactant_indices_val = [1]
product_indices_val = [3, 4]
reactants_stoichs_val = [1]
products_stoichs_val = [1, 2]

n_species_val = length(stoich_coeffs_val)
n_reactants_val = 1
n_products_val = 2

rows = 4
n_nodes = rows

heat_capacity_vec_val = [20, 20, 20, 20, 20]
reactants_stoichs_val = [1]
products_stoichs_val = [1, 2]

@named left_bcs = SpeciesSource(n_species=n_species_val, n_reactants=n_reactants_val, n_products=n_products_val, F_fixed=[1.0, 1.0, 0.01, 0.01, 0.01], T_fixed=500.0, p_fixed=5e5)

ReactorSegment(name=Symbol("seg"),
                n_species = n_species_val,
                n_reactants = length(reactants_stoichs_val),
                n_products = length(products_stoichs_val),
                catalyst_weight = 5,
                heat_of_reaction = 90.7e3, # J/mol
                heat_capacity_vec = [81.6, 33.6, 29.1, 28.8, 37.1], # J/(mol*K)
                stoich_coeffs = stoich_coeffs_val,
                kinetics = ReactionKinetics(name=Symbol("kinetics"), 
                                        n_species = n_species_val,
                                        n_reactants = length(reactants_stoichs_val),
                                        n_products = length(products_stoichs_val),
                                        #reactant_indices=reactant_indices_val,
                                        #product_indices=product_indices_val,
                                        reactants_stoichs=reactants_stoichs_val,
                                        products_stoichs=products_stoichs_val,
                                        kf_A=1.0e5, kf_Ea=8.0e4, kr_A=1.0e12, kr_Ea=1.5e5)
)

segments = [ReactorSegment(name=Symbol("seg", i),
                   n_species=n_species_val,
                   n_reactants=length(reactant_indices_val),
                   n_products=length(product_indices_val),
                   catalyst_weight = catalyst_per_segment,
                   heat_of_reaction = 90.7e3, # J/mol
                   heat_capacity_vec = [81.6, 33.6, 29.1, 28.8, 37.1], # J/(mol*K)
                   stoich_coeffs = stoich_coeffs_val,
                   kinetics = ReactionKinetics(name=Symbol("kinetics", i), reactants_stoichs=reactants_stoichs_val,
                                              products_stoichs=products_stoichs_val,
                                              reactant_indices=reactant_indices_val,
                                              product_indices=product_indices_val,
                                              kf_A=1.0e5, kf_Ea=8.0e4, kr_A=1.0e12, kr_Ea=1.5e5)
                   ) for i in 1:rows]

#nodes = [ReactorSegment(name=Symbol("Reactor_Segment_", i), catalyst_weight=1, heat_of_reaction=-92.0, heat_capacity_vec=heat_capacity_vec_val, reactants_stoichs=reactants_stoichs_val, products_stoichs_vec=products_stoichs_vec_val) for i in 1:rows]


#!!!!IMPORTANT!!!! NOTE: For some reason, whenever you do F(t)[1:n] YOU CANNOT DO A DESCRIPTION !!!!IMPORTANT!!!!


right_bcs = SpeciesSource(F_fixed=[1.0, 1.3, 0.001, 0.001, 0.001], T_fixed=293.15, p_fixed=101325, name=Symbol("Species_Source"))

connections = Equation[]

push!(connections, connect(left_bcs.port, nodes[1].port))
[push!(connections, connect(nodes[i].port, nodes[i].port)) for i in 1:(rows-1)]
push!(connections, connect(nodes[end].port, right_bcs.port))

# --- Simulation Setup ---
eqs_total = System[]
all_systems = vcat(
    vec(nodes),
    vec(left_bcs),
    vec(right_bcs)
)
all_systems

@named rod = ODESystem(connections, t, systems=all_systems)

sys = structural_simplify(rod)  

tspan = (0.0, 10.0)

prob = ODEProblem(sys, [], tspan)

sol = solve(prob)
