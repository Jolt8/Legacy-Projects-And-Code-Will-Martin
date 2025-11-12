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
    end
    
    @variables begin
        F(t)[1:n_species]#, [description="test"] #description breaks this for some reason
        T(t)
        p(t)
        rate(t)
        #forward_term(t)
        #reverse_term(t)
    end

    @parameters begin
        reactants_stoichs[1:n_species]
        products_stoichs[1:n_species]
        kf_A, [description="Forward Arrhenius pre-exponential factor"]
        kf_Ea, [description="Forward reaction activation energy"]
        kr_A, [description="Reverse Arrhenius pre-exponential factor"]
        kr_Ea, [description="Reverse reaction activation energy"]
        R_gas = 8.314
        
        F_total = sum(F)
        
        kf = kf_A * exp(-kf_Ea / (R_gas * T))
        kr = kr_A * exp(-kr_Ea / (R_gas * T))
        
        partial_pressures[1:n_species] = [p * F[i] / F_total for i in 1:n_species] #this is the only syntax that works for doing arrays within parameters
        
        pre_f_t[1:n_species] = [partial_pressures[i]^reactants_stoichs[i] for i in 1:n_species] #this is a pretty brittle way because we're relying on num^0 = 1, Too Bad!
        pre_r_t[1:n_species] = [partial_pressures[i]^products_stoichs[i] for i in 1:n_species]
        forward_term = prod(pre_f_t)
        reverse_term = prod(pre_r_t)
        #NOTE: WHAT THE FUCK: Why does the forward and reverse only work if we do a "pre" version of them
    end
    #NOTE: HOLD UP, if you define the variables before the parameters, you can do sum(F)
    
    @equations  begin
        rate ~ kf * forward_term - kr * reverse_term
    end
end
"""@named test = ReactionKinetics(name=Symbol("kinetics", 1), 
                                n_species=5,
                                n_reactants=1,
                                n_products=2,
                                reactants_stoichs=[1, 0, 0, 0, 0],
                                products_stoichs=[0, 0, 1, 2, 0],
                                kf_A=1.0e5, kf_Ea=8.0e4, kr_A=1.0e12, kr_Ea=1.5e5
)
unknowns(test)
parameters(test)
test.pre_f_t[1]"""

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
        D(p) ~ -1 * ((150 * μ * G*(1 - ϵ)^2) / (Dp^2 * ϵ^3) + (1.75 * G^2 * (1 - ϵ)) / (Dp * ϵ^3))
    end
end

#Note sure if a reactor segment requires a 1D Element to connect them 
#Perhaps this would be a good use case for trying out the finite volume method in mtkcomponents
@mtkmodel ReactorSegment begin
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
        
        reactants_stoichs[1:n_species]
        products_stoichs[1:n_species]
        kf_A
        kf_Ea
        kr_A
        kr_Ea
    end
    
    @components begin
        inlet = ReactorPort(n_species=n_species, n_reactants=n_reactants, n_products=n_products)
        outlet = ReactorPort(n_species=n_species, n_reactants=n_reactants, n_products=n_products)
        kinetics = ReactionKinetics(
            n_species=n_species, 
            n_reactants=n_reactants, 
            n_products=n_products,
            reactants_stoichs=reactants_stoichs,
            products_stoichs=products_stoichs,
            kf_A=kf_A, 
            kf_Ea=kf_Ea, 
            kr_A=kr_A, 
            kr_Ea=kr_Ea)
        #pressure_drop = ErgunPressureDrop()
    end


    @variables begin
        F(t)[1:n_species]
        T(t) = 293.13
        p(t) = 101325
    end

    @equations begin
        kinetics.F ~ F
        kinetics.T ~ T
        kinetics.p ~ p
        #pressure_drop.p ~ p

        D(F) ~ [kinetics.rate * stoich_coeff for stoich_coeff in stoich_coeffs] #don't use dot here: [kinetics.rate .* stoich_coeffs]

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
reactants_stoichs_val = [1, 0, 0, 0, 0]
products_stoichs_val = [0, 0, 1, 2, 0]

@named left_bcs = SpeciesSource(n_species=n_species_val, n_reactants=n_reactants_val, n_products=n_products_val, F_fixed=[1.0, 1.0, 0.01, 0.01, 0.01], T_fixed=500.0, p_fixed=5e5)

#NOTE: whenever it says that it cannot convert a type into an equation, it usually means you are incorrectly defining equation as not an equation somewhere in your model
ReactorSegment(name=Symbol("seg"),
    n_species = n_species_val,
    n_reactants = 1,
    n_products = 1,
    catalyst_weight = 5,
    heat_of_reaction = 90.7e3, # J/mol
    heat_capacity_vec = [81.6, 33.6, 29.1, 28.8, 37.1], # J/(mol*K)
    stoich_coeffs = stoich_coeffs_val,
    reactants_stoichs=reactants_stoichs_val,
    products_stoichs=products_stoichs_val,
    kf_A=1.0e5, kf_Ea=8.0e4, kr_A=1.0e12, kr_Ea=1.5e5)


segments = [ReactorSegment(name=Symbol("seg", i),
    n_species = n_species_val,
    n_reactants = 1,
    n_products = 1,
    catalyst_weight = 5,
    heat_of_reaction = 90.7e3, # J/mol
    heat_capacity_vec = [81.6, 33.6, 29.1, 28.8, 37.1], # J/(mol*K)
    stoich_coeffs = stoich_coeffs_val,
    reactants_stoichs=reactants_stoichs_val,
    products_stoichs=products_stoichs_val,
    kf_A=1.0e5, kf_Ea=8.0e4, kr_A=1.0e12, kr_Ea=1.5e5) for i in 1:rows
]

#NOTE: in the future, if you want to change a value like kf_A halfway through the reactor just change kinetics directly like: seg.kinetics.kf_A

#nodes = [ReactorSegment(name=Symbol("Reactor_Segment_", i), catalyst_weight=1, heat_of_reaction=-92.0, heat_capacity_vec=heat_capacity_vec_val, reactants_stoichs=reactants_stoichs_val, products_stoichs_vec=products_stoichs_vec_val) for i in 1:rows]

#

#!!!!IMPORTANT!!!! NOTE: For some reason, whenever you do F(t)[1:n] YOU CANNOT DO A connect = Flow !!!!IMPORTANT!!!!
#actually, I'm not sure if the above is true, this requires more testing in a later complete system

@named right_bcs = SpeciesSource(n_species=n_species_val, n_reactants=n_reactants_val, n_products=n_products_val, F_fixed=[1.0, 1.0, 0.01, 0.01, 0.01], T_fixed=500.0, p_fixed=5e5)

connections = Equation[]

push!(connections, connect(left_bcs.port, segments[1].inlet))
[push!(connections, connect(segments[i].outlet, segments[i+1].inlet)) for i in 1:(rows-1)]
push!(connections, connect(segments[end].outlet, right_bcs.port))

# --- Simulation Setup ---
eqs_total = System[]
all_systems = vcat(
    segments,
    left_bcs,
    right_bcs
)
all_systems

@named rod = ODESystem(connections, t, systems=all_systems)

sys = structural_simplify(rod)  

tspan = (0.0, 10.0)

prob = ODEProblem(sys, [], tspan)

sol = solve(prob)
