# Model Borders (or: On Voltage and Current Sources)

In this document we'll go through some very important concepts relating to the intricate interplay between ModelingToolkit and NetworkDyanmics.jl.
Make sure to read the PowerDynamics docs on [Modeling Concepts](@ref) first.
Also check out the NetworkDyanmics.jl docs on the [Mathematical Model](@extref), which shows in detail what will be recapped here in short form.

Since handling huge symbolic models can be computationaly intensive, the whole point of this library to have a clear distinction: truly acausal symbolic models for the indivdual components (like Generators, Loads, RES, ...) which gets then symbolicially simplified and compiled into [`VertexModel`](@extref NetworkDyanmics.VertexModel) and [`EdgeModel`](@extref NetworkDyanmics.EdgeModel) objects, which are the actual models used for the network simulation.
Those models have a clear input-output structure related to the interconnection of potentials and flows on the network.
Namely:
- VertexModels sit at the buses. As an input, they see the current sum of all connected lines. Their job is to establish a voltage at that point.
- EdgeModels sit on the edges. As an input they see the voltages on both ends. Their job is to establish the current flow on that line.

```asciiart
                                 more edges
                                     △
n ⋯───╮             ╭────────────────┼────────────────╮             ╭───⋯ n
e     │             │        votlage │ u out          │             │     e
x  ┏━━▽━━━━━━━━━━━━━▽━━┓   ╔═════════△═════════╗   ┏━━▽━━━━━━━━━━━━━▽━━┓  x
t  ┃ EdgeModel         ┃   ║ VertexModel       ║   ┃ EdgeModel         ┃  t
   ┃ ẋ = f(x, u, p, t) ┃   ║ ẋ = f(x, i, p, t) ║   ┃ ẋ = f(x, u, p, t) ┃
n  ┃ i = g(x, u, p, t) ┃   ║ u = g(x, p, t)    ║   ┃ Φ = g(x, u, p, t) ┃  n
o  ┗━━▽━━━━━━━━━━━━━▽━━┛   ╚═════════△═════════╝   ┗━━▽━━━━━━━━━━━━━▽━━┛  o
d     │     current │ i out        ╭─┴─╮      current │ i out       │     d
e ⋯───╯             ╰──────────────▷ + ◁──────────────╯             ╰───⋯ e
                                   ╰─△─╯
                                     │
                                 more edges
```

So at the core, a VertexModel **musst act like a voltage source** (we'll come to the excption soon), while EdgeModels **musst act like a current source**.

You might ask: what happens if that is not the case? For example, at typical generator model
might be modeled as some kind of voltage source behind a stator resistance/reactance.


```asciiart
       u_device    Susceptance  u_bus
                ●──────███───→───●──←───
                │   i_device        i_grid
voltage source (↗)
                │
                │
                ⏚
```
The obvious equations for this system are
```math
\begin{algined}
   i_{\mathrm{device}} &= Y\,u_{\mathrm{bus}}
   i_{\mathrm{grid}} &= i_{\mathrm{device}}
\end{aligned}
```
which looks an aweful lot like a system with voltage as the natural input and current as the output.
However, we can use Kirchoff-Current-Law (KCL) to resolve this problem, by introducing a
constraint equation:
```math
\begin{algined}
   i_{\mathrm{device}} &= Y\,u_{\mathrm{bus}}
   0  &= i_{\mathrm{device}} - i_{\mathrm{grid}}
\end{aligned}
```
Here, even though the second equation looks similar it is not a straight forward relation to calculated $i_{\mathrm{device}} as an output. Instead it is an **fully implicit equation** for bus voltage $u_{\mathrm{bus}}$. I.e. in closed loop with the Network, the solver has to find a $u_{\mathrm{bus}}$ such that the KCL is fulfilled.
This is kind of modeling is quite typical for power grid simulations.

However, we have two practical problem arising from merging all devices connected to a bus into a single VertexModel and using KCL as an output constraint:
1. Large Models: In some simulations you might encounter large aggregated models, like ten generators including controller dynamics at a single bus. This can lead to single bus models with hundreds of states, leading to long model compilation times.
2. EMT Models: If you happen to do EMT modeling you currents mostly become differential states. Those **must** not be coupled using KCL.

Problem 1 can be solved using special [Current Injector Bus](@ref) models described below.
Problem 2 is slighly more fundamental and will be adressed in detail first.

## Challenges of EMT Models 
In EMT style modeling everything looks a bit different. All of the sudden you get differentail equations for voltages and currents. Lets consider the voltage-source-behind-susceptance model again:
```math
\begin{algined}
   \frac{\mathrm{d}i_{\mathrm{device}}}{\mathrm{d}t} &= f_Y(u_{\mathrm{bus}})
   0  &= i_{\mathrm{device}} - i_{\mathrm{grid}}
\end{aligned}
```
Lets also assume that you line is something like an EMT RL model now, which also defines its current as a differential state.
Now solving the KCL becomes a real problem!
One way to intuitively think about DAEs is, that you want to locally "twiddle" the algebraic states (in this case our bus voltage) until you solve the constraint equations while the differential states remain fixed. With algebraic current equations, that is possible.
Once both $i_{\mathrm{device}}$ and $i_{\mathrm{grid}}$ are differential states, that is no longer possible. There is no way to find $u_{\mathrm{bus}}$ such that the KCL is fulfilled.
From a mathematical standpoint, we created an higher index DAE.

There are two ways to solve this problem:
1. Replace the *algebraic* state by a *differential* state
2. Introduce some *algebraic* current to the KCL.

The first solution is to consider some capacity at the bus:
```math
\begin{algined}
   \frac{\mathrm{d}i_{\mathrm{device}}}{\mathrm{d}t} &= f_Y(u_{\mathrm{bus}})
   C\,\frac{\mathrm{d}u_{\mathrm{bus}}} &= i_{\mathrm{device}} - i_{\mathrm{grid}}
\end{aligned}
```
Now, $u_{\mathrm{bus}}$ is a differential state and there is no KCL constraint anymore.
Explicitly modeling this capacitance is close to the true physical reality of electrical systems. 
In some sense the constraint only arived after we went to the limit of $C \to 0$ after discretizing the [Telegrapher's Equation](https://en.wikipedia.org/wiki/Telegrapher%27s_equations) of the actual conductor.

The second solution is to introduce some algebraic current into the KCL:
```math
\begin{algined}
   \frac{\mathrm{d}i_{\mathrm{device}}}{\mathrm{d}t} &= f_Y(u_{\mathrm{bus}})
   0 &= i_{\mathrm{device}} - i_{\mathrm{grid}} + f_i(u_{\mathrm{bus}})
\end{aligned}
```
This can be achieved by connecting a resistor to ground for example.
Alternatively, this is implicily achieved if $i_{\mathrm{grid}}$ is not a pure differential state, but rather defined by some algebraic equation.
This is trivially fulfilled for standard piline phasor models. Alternatviely one could add static R shunts to the line models.
If this is the case $i_{\mathrm{grid}} = f_i(u_{\mathrm{bus}})$ itself directly depens on the busbar voltage, making the KCL solvable again. 
In some sense you could say, that the edge model provides direct FF coupling of vertex output $u_{\mathrm{bus}}$ to vertex intput $i_{\mathrm{grid}}$.

!!! tip: Tell NetworkDynamics to assume static lines
    Sometimes you get problems on `compile_bus` or `VertexModel` with error messages like "to many highest order equations" or errors related to required input derivatives.
    This happens if try to compile a bus model with dynamic state $i_{\mathrm{device}}$. You can use the `assume_io_coupling=true` keyword to make MTK aware of the direct feedback relation
$i_{\mathrm{grid}} = f_i(u_{\mathrm{bus}})$. 

We explained the problem at the example of the voltage soure behind susceptance model.
A similar problem arises if you try to define a dynamic PiLine for example. You can't do that, because then you would define separate differential equations for a single state (two connected PI lines would define differential equations for the same bus voltage). Circumvent this by defining dynamic Tau-lines or aggregating the dynamic shunts at the vertex models.
    
Sometimes, it is not nice to alter your natural vertex current sources by pulling helper elements like shunts into the vertex model equations. For that usecase we have a special instrument: Curreent Injector Buses.


## Current Injector Bus 
Sometimes you just want to hook multiple current injectors on a single busbar.
NetworkDynamics allows you to do so by introducing a special way of modeling "vertex clusters".
A vertex cluster contains of a single hub vertex which acts like a normal voltage source element.
However you are allowd to attach multiple other current-source-like vertices to the hub using
so called [`LoobpackConnection`](@extref NetworkDynamics.LoobpackConnection) edges.

A loopback conenction is a special edge which **directly** connects vertices (similar to a zero impedance line).
This means the input-output system of the satelites is reversed: they see the hub voltage as input and inject current into the hub bus.

```asciiart
                  Hub    Loopback  Satelites
                ╭──────╮╭────────╮╭──────────╮
      
  ┏━━━━━━┓                        ┏━━━━━━━━━━┓
⋯─┨Line A┠──╮   ┏━━━━━━┓    ╭─────┨Injector A┃
  ┗━━━━━━┛  ╰───┨      ┠────╯     ┗━━━━━━━━━━┛
  ┏━━━━━━┓      ┃ Σi=0 ┃
⋯─┨Line B┠──────┨      ┠────╮     ┏━━━━━━━━━━┓
  ┗━━━━━━┛      ┗━━━━━━┛    ╰─────┨Injector B┃
                                  ┗━━━━━━━━━━┛

                ╰────────────────────────────╯
                     Vertex Cluster
```

Often, the Hub bus will be a pure junction bus (pure KCL constraint).
For EMT sims it might be a pure capacitive shunt as described above.
In theory, it can be any type of model what soever.


### Current Injector Bus Example
We want to model a simple system with two generators connected to a load via a PiLine.

We'll model this system in two different ways: 
First we model it by merging two generators into a single vertex model using KCL on the interconnection.
Then we'll model the system as a 4-Bus model using loopback connections.

```asciiart
       1
G1 (~)─┨          2
       ┠──────────╂─▷ L
G2 (~)─┨    
```

```asciiart
       1   3
G1 (~)─╊═══┫      4
           ╂──────╂─▷ L
G2 (~)─╊═══┫
       2    
```

#### Modeling as normal 2-Bus System
```@example
using ModelingToolkit, PowerDynamics, NetworkDynamics, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve, CairoMakie
@named load = ConstantYLoad()
@named loadbus = compile_bus(MTKBus(load); pf=pfPQ(P=-1, Q=-0.1))
```
Lets define a quick perturbance event by changing the load power at t=1s:
```@example
affect = ComponentAffect([], [:load₊G]) do u, p, ctx
    if ctx.t == 0.1
       p[:load₊G] *= 100 # admittance increase by 10%
    elseif ctx.t == 0.2
       p[:load₊G] /= 100 # admittance increase by 10%
    end
end
cb = PresetTimeComponentCallback([0.1, 0.2], affect)
add_callback!(loadbus, cb)
nothing #hide
```

Also lets define the line model
```@example
@named line = compile_line(MTKLine(PiLine(X=0.1, R=0.01;name=:pi)))
```

Next we define the mixed generaor bus:
```@example
genp = (vf_input=false, τ_m_input=false, S_b=100, V_b=1, ω_b=2π*50, R_s=0.000125, T″_d0=0.01, T″_q0=0.01, X_ls=0.01460, X_d=0.1460, X′_d=0.0608, X″_d=0.06, X_q=0.1000, X′_q=0.0969, X″_q=0.06, T′_d0=8.96, T′_q0=0.310, H=23.64)

@named genA = SauerPaiMachine(; genp...)
@named genB = SauerPaiMachine(; genp...)
genABmod = MTKBus(genA, genB; name=:GEN1)
@named genABbus = compile_bus(genABmod, pf=pfSlack(V=1))
```
Here we actually hit a problem the first time.
Our combined bus has the single powerflow model of a slack.
However during initialization, we know nothing about the powersharing!
What we want, is that the first generator acts as a slack and the second generator acts as a
PV node.
We can achieve this behavior by defining additional init constraints:

```@example
psharing = @initconstraint begin
    :genB₊P - 0.45 # V = 1
    :genB₊Q - 0.1 # V = 1
end
set_initconstraint!(genABbus, psharing)
```

then we can build the network
```@example
line_ab = EdgeModel(line; src=:genABbus, dst=:loadbus)
nw = Network([genABbus, loadbus], [line_ab])
s0 = initialize_from_pf(nw)
prob = ODEProblem(nw, s0, (0.0, 3.0))
sol = solve(prob, Rodas5P());
let
    fig = Figure()
    ax1 = Axis(fig[1,1], title="Generator Voltage")
    lines!(ax1, sol, idxs=VIndex(:genABbus, :genA₊v_mag), color=Cycled(1), label="Gen A")
    lines!(ax1, sol, idxs=VIndex(:genABbus, :genB₊v_mag), color=Cycled(2), label="Gen B")
    axislegend(ax1)
    ax2 = Axis(fig[2,1], title="Generator Active Power P")
    lines!(ax2, sol, idxs=VIndex(:genABbus, :genA₊P), color=Cycled(1))
    lines!(ax2, sol, idxs=VIndex(:genABbus, :genB₊P), color=Cycled(2))
    ax3 = Axis(fig[3,1], title="Generator Reactive Power Q")
    lines!(ax3, sol, idxs=VIndex(:genABbus, :genA₊Q), color=Cycled(1))
    lines!(ax3, sol, idxs=VIndex(:genABbus, :genB₊Q), color=Cycled(2))
    fig
end
```

No we model the same system using loopback


```@example
@named gen = SauerPaiMachine(; genp...)
@named genAbus = compile_bus(MTKBus(gen), current_source=true)
```
note how we used current source to generate a model with u input i output
make sure that you also compile the powerflow model as a current source
```@example
set_pfmodel!(genAbus, pfSlack(V=1; current_source=true))
# claude: introuduce the trick here from EMT Toyexample where we show the error something with try cath and # hide and maybe showerror
```
Oh no that failed! The reason is similar to the one described above. The slack is modeld as a constraint `u = u_set`. This cannot be resolved unless we assume instantaneous feedback from bus voltage to bus current. We can use the `assume_io_coupling=true` trick to fix this:
```@example
set_pfmodel!(genAbus, pfSlack(V=1; current_source=true, assume_io_coupling=true))
nothing # hide
```

similar, we do it fo gen B but this time with a PQ model
```@example
@named genBbus = compile_bus(MTKBus(gen), current_source=true)
set_pfmodel!(genBbus, pfPQ(P=0.45, Q=0.1, current_source=true))
nothing #hide
```
lastly we need to define the junction bus
```@example
@named junction = compile_bus(MTKBus()) # pur KCL
```
to connect we define our two loopback connections
```@example
loopbackA = LoopbackConnection(; src=:genAbus, dst=:junction, potential=[:u_r, :u_i], flow=[:i_r, :i_i])
loopbackB = LoopbackConnection(; src=:genBbus, dst=:junction, potential=[:u_r, :u_i], flow=[:i_r, :i_i])
```

```@example
line_junction = EdgeModel(line; src=:junction, dst=:loadbus)
nw = Network([genAbus, genBbus, junction, loadbus], [loopbackA, loopbackB, line_junction])
```
unsurprisingly now we have a network with 4 vertices and 3 edges
```@example
s0 = initialize_from_pf(nw; tol=1e-9, nwtol=1e-8)
```
The network still initialies fine. Since we have separate powerflow models for each generator now the powersharing is directly defined by the powerflow models.

```@example
prob = ODEProblem(nw, s0, (0.0, 3.0))
sol = solve(prob, Rodas5P());
let
    fig = Figure()
    ax1 = Axis(fig[1,1], title="Generator Voltage")
    lines!(ax1, sol, idxs=VIndex(:genAbus, :gen₊v_mag), color=Cycled(1), label="Gen A")
    lines!(ax1, sol, idxs=VIndex(:genBbus, :gen₊v_mag), color=Cycled(2), label="Gen B")
    axislegend(ax1)
    ax2 = Axis(fig[2,1], title="Generator Active Power P")
    lines!(ax2, sol, idxs=VIndex(:genAbus, :gen₊P), color=Cycled(1))
    lines!(ax2, sol, idxs=VIndex(:genBbus, :gen₊P), color=Cycled(2))
    ax3 = Axis(fig[3,1], title="Generator Reactive Power Q")
    lines!(ax3, sol, idxs=VIndex(:genAbus, :gen₊Q), color=Cycled(1))
    lines!(ax3, sol, idxs=VIndex(:genBbus, :gen₊Q), color=Cycled(2))
    fig
end
```
