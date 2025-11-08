using ModelingToolkit
using ModelingToolkitStandardLibrary.Thermal: HeatPort, Element1D
using DifferentialEquations
using Plots

import ModelingToolkit: t_nounits as t, D_nounits as D

# A simple 0D heat capacitor, representing one node.
@mtkmodel HeatCapacitor begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        C = 1.0, [description = "Heat capacity of element"]
    end
    @variables begin
        T(t) = 293.15, [description = "Temperature of element [K]"]
        der_T(t) = 0.0, [description = "Time derivative of temperature [K/s]"]
    end

    @equations begin
        T ~ port.T
        D(T) ~ port.Q_flow / C
    end
end

# A 0D thermal conductor, representing the connection between nodes.
@mtkmodel ThermalConductor begin
    @extend Q_flow, dT = element1d = Element1D()
    @parameters begin
        G = 1.0, [description = "Thermal conductance [W/K]"]
    end
    @equations begin
        Q_flow ~ G * dT
    end
end

# A constant temperature boundary condition
@mtkmodel SomeFixedTemperature begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        T_fixed = 303.15, [description = "Fixed temperature [K]"] # 30 °C
    end
    @equations begin
        port.T ~ T_fixed
    end
end

# A constant heat flow boundary condition
@mtkmodel SomeHeatSource begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        Q_flow_fixed = 10.0, [description = "Fixed heat flow [W]"]
    end
    @equations begin
        port.Q_flow ~ Q_flow_fixed
    end
end

# Now, build the 1D PDE system (the rod) by connecting them!
n = 3

# Instantiate all the components
@named left_boundary = SomeHeatSource(Q_flow_fixed = 5.0)
@named right_boundary = SomeFixedTemperature(T_fixed = 293.15)
nodes = [HeatCapacitor(name=Symbol("node", i), C=1/n) for i in 1:n]
conductors = [ThermalConductor(name=Symbol("cond", i), G=n) for i in 1:(n-1)]

# 2. Create connections programmatically
connections = Equation[]
# Connect internal nodes and conductors
for i in 1:(n-1)
    push!(connections, connect(nodes[i].port, conductors[i].port_a))
    push!(connections, connect(conductors[i].port_b, nodes[i+1].port))
end
# Connect boundaries
push!(connections, connect(left_boundary.port, nodes[1].port))
push!(connections, connect(right_boundary.port, nodes[n].port))

# 3. Create the ODESystem
all_systems = vcat(nodes, conductors, left_boundary, right_boundary)
@named rod = ODESystem(connections, t, systems=all_systems)

#We should probably get the function
"""
@mtkmodel HeatRod begin
    @components begin
        # Create n capacitors (the nodes)
        nodes = [HeatCapacitor(C=1/n) for i in 1:n]
        # Create n-1 conductors (the connections)
        conductors = [Element1D() for i in 1:(n-1)]

        # Add boundary condition components
        left_boundary = SomeHeatSource()
        right_boundary = SomeFixedTemperature(T_fixed=293.15) # 20 °C
    end
    
    @equations begin
        # This loop IS the discretization of the PDE
        for i in 1:(n-1)
            connect(nodes[i].port, conductors[i].port_a)
            connect(conductors[i].port_b, nodes[i+1].port)
        end
        
        # Connect the ends of the rod to the boundary conditions
        connect(left_boundary.port, nodes[1].port)
        connect(right_boundary.port, nodes[n].port)
    end
end
"""
# --- Simulation Setup ---

# Instantiate the model
#@named rod = HeatRod()

# Create a system and simplify the equations
sys = structural_simplify(rod)

# Define the simulation time span
tspan = (0.0, 10.0)

# Create an ODEProblem
prob = ODEProblem(sys, [], tspan)

# Solve the problem
sol = solve(prob)

# --- Plotting Results ---

# Extract the temperature variables for each node
#temp_vars = [rod.nodes[i].T for i in 1:10]

# Plot the solution
plot(sol,
     vars=[node.T for node in nodes],
     title="Temperature Distribution in Heat Rod",
     xlabel="Time (s)",
     ylabel="Temperature (K)",
     legend=:topleft,
     #labels=["Node " .* string.(1:10) for i in 1:10]'
    )