using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t, D_nounits as D

# Basic electric components
@connector function Pin(; name)
    @variables v(t)=1.0 i(t)=1.0 [connect=Flow]
    System(Equation[], t, [v, i], [], name = name)
end

function Ground(; name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(System(eqs, t, [], [], name = name), g)
end

function ConstantVoltage(; name, V = 1.0)
    val = V
    @named p = Pin()
    @named n = Pin()
    @parameters V = V
    eqs = [V ~ p.v - n.v
           0 ~ p.i + n.i]
    compose(System(eqs, t, [], [V], name = name), p, n)
end

@connector function HeatPort(; name)
    @variables T(t)=293.15 Q_flow(t)=0.0 [connect=Flow]
    System(Equation[], t, [T, Q_flow], [], name = name)
end

function HeatingResistor(; name, R = 1.0, TAmbient = 293.15, alpha = 1.0)
    @named p = Pin()
    @named n = Pin()
    @named h = HeatPort()
    @variables v(t) RTherm(t)
    @parameters R=R TAmbient=TAmbient alpha=alpha
    eqs = [RTherm ~ R * (1 + alpha * (h.T - TAmbient))
           v ~ p.i * RTherm
           h.Q_flow ~ -v * p.i # -LossPower
           v ~ p.v - n.v
           0 ~ p.i + n.i]
    compose(System(eqs, t, [v, RTherm], [R, TAmbient, alpha],
            name = name), p, n, h)
end

function HeatCapacitor(; name, rho = 8050, V = 1, cp = 460, TAmbient = 293.15)
    @parameters rho=rho V=V cp=cp
    C = rho * V * cp
    @named h = HeatPort()
    eqs = [
        D(h.T) ~ h.Q_flow / C
    ]
    compose(System(eqs, t, [], [rho, V, cp],
            name = name), h)
end

function Capacitor(; name, C = 1.0)
    @named p = Pin()
    @named n = Pin()
    @variables v(t) = 0.0
    @parameters C = C
    eqs = [v ~ p.v - n.v
           0 ~ p.i + n.i
           D(v) ~ p.i / C]
    compose(System(eqs, t, [v], [C],
            name = name), p, n)
end

function parallel_rc_model(i; name, source, ground, R, C)
    resistor = HeatingResistor(name = Symbol(:resistor, i), R = R)
    capacitor = Capacitor(name = Symbol(:capacitor, i), C = C)
    heat_capacitor = HeatCapacitor(name = Symbol(:heat_capacitor, i))

    rc_eqs = [connect(source.p, resistor.p)
              connect(resistor.n, capacitor.p)
              connect(capacitor.n, source.n, ground.g)
              connect(resistor.h, heat_capacitor.h)]

    compose(System(rc_eqs, t, name = Symbol(name, i)),
        [resistor, capacitor, source, ground, heat_capacitor])
end

V = 2.0
@named source = ConstantVoltage(V = V)
@named ground = Ground()
N = 50
Rs = 10 .^ range(0, stop = -4, length = N)
Cs = 10 .^ range(-3, stop = 0, length = N)
rc_systems = map(1:N) do i
    parallel_rc_model(i; name = :rc, source = source, ground = ground, R = Rs[i], C = Cs[i])
end;
@variables E(t) = 0.0
eqs = [
    D(E) ~ sum(((i, sys),) -> getproperty(sys, Symbol(:resistor, i)).h.Q_flow,
    enumerate(rc_systems))
]
@named _big_rc = System(eqs, t, [E], [])
@named big_rc = compose(_big_rc, rc_systems)

sys = structural_simplify(big_rc)

tspan = (0.0, 10.0)

prob = ODEProblem(sys, [], tspan)

sol = solve(prob) #bro even this fails