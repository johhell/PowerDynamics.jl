#=
# [Linear Analysis of a 4-Bus System](@id linear-analysis)

This example can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

This example demonstrates **eigenvalue analysis** and **impedance-based Bode analysis** of a
4-bus EMT power system.

!!! tip "Replicated SimplusGT Example"
    The system parameters and topology are taken from [SimplusGT](https://github.com/Future-Power-Networks/Simplus-Grid-Tool) (Simplus Grid Tool).
    To compare the results this example reimplements their basic 4-bus default example.
    Since we use the same models, same parameters and so on we can directly compare our results to the SimplusGT results.

    Please note that this tutorial only replicates a small part of what SimplusGT is capable of,
    so please check out their great toolbox for in-depth linear analysis of powersystems!
=#


#=
## Synchronous Machine with Stator Dynamics
In order to replicate the results from SimplusGT we need to use identical or nearly identical models to
SimplusGT. Most of the models needed are present in our Library, we still need to implement a custom machine model though.
Since modeling is not the focus on this tutorial, it is provide without any further explaination.

```@raw html
<details class="admonition is-details">
<summary class="admonition-header">Definition of custom Machine Model</summary>
<div class="admonition-body">
SimplusGT uses a special machine model for their Type 0 Apparatus which:
- retains the stator flux dynamics (i.e. defines the current output as an differential equation),
- assumes constant field flux (parameter which is initialized at operation point) and
- contains no governor dynamics (i.e. mechanical torque is a constant parameter initialized at operation point).
```
=#
using PowerDynamics
using PowerDynamics.Library
using PowerDynamics.Library.ComposableInverter
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using NetworkDynamics
using CairoMakie

@mtkmodel SyncMachineStatorDynamics begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        J, [description="Inertia constant [MWs²/MVA]"]
        D, [description="Damping coefficient [pu]"]
        wL, [description="Stator inductance * base frequency [pu]"]
        R, [description="Stator resistance [pu]"]
        ω0, [description="Base frequency [rad/s]"]
        psi_f, [guess=1, description="Field flux linkage [pu]"]
        T_m, [guess=1, description="Mechanical torque [pu]"]
    end
    @variables begin
        i_d(t), [guess=0, description="d-axis stator current [pu]"]
        i_q(t), [guess=1, description="q-axis stator current [pu]"]
        w(t), [guess=2*pi*50, description="Rotor speed deviation [rad/s]"]
        theta(t), [guess=0, description="Rotor angle [rad]"]
        v_d(t), [guess=1, description="d-axis terminal voltage [pu]"]
        v_q(t), [guess=0, description="q-axis terminal voltage [pu]"]
        psi_d(t), [guess=1, description="d-axis flux linkage [pu]"]
        psi_q(t), [guess=0, description="q-axis flux linkage [pu]"]
        Te(t), [guess=1, description="Electrical torque [pu]"]
    end
    begin
        J_pu = J*2/ω0^2
        D_pu = D/ω0^2
        L = wL/ω0
        T_to_loc(α)  = [ cos(α) sin(α);
                        -sin(α)  cos(α)]
        T_to_glob(α) = T_to_loc(-α)
    end
    @equations begin
        # output transformation (global dq/local dq)
        [terminal.i_r, terminal.i_i] .~ -T_to_glob(theta)*[i_d, i_q]
        [v_d, v_q] .~ T_to_loc(theta)*[terminal.u_r, terminal.u_i]
        # electromechical equations
        psi_d ~ L*i_d
        psi_q ~ L*i_q - psi_f
        Te ~ psi_f * i_d
        Dt(i_d) ~ (v_d - R*i_d + w*psi_q)/L
        Dt(i_q) ~ (v_q - R*i_q - w*psi_d)/L
        # swing equation
        Dt(w) ~ (Te - T_m - D_pu*w)/J_pu
        Dt(theta) ~ w - ω0
    end
end
nothing #hide #md

#=
```@raw html
</div>
</details>
```

## 4-Bus Network Setup
The topolgoy of the network looks like this:
```asciiart
       (~)
    2╺┯━┷━┯╸
      │   │   (GFL)
1╺┯━━━┿╸  │   ╺━┿━╸4
 (~)  │   │╭────╯
     ╺┷━┯━┷┷╸3
      (GFM)

```
We have 4 Buses in total, two of which have a synchronous machine connected, one has
a grid forming droop inverter and one has a grid following inverter.

!!! note "EMT Models"
    The entire modeling in SimplusGT uses **EMT Components**, meaning that the powerlines
    are modeled using [`DynamicRLBranch`](@ref) models. At each of the buses there is
    shunt consisting of a parallel capacitor and resistor.

Since the devices are all natural current sources and we have explicit shunts at all
buses (self-edges in SimplusGT data), we use current source modeling accoding to the docs [on voltage and current sources](@ref vc-and-cs)
Specificially, we model each bus using the RC-shunt.
The devices are connected to those "shut buses" as current sources via [`LoopbackConnection`](@extref NetworkDynamics.LoopbackConnection) components.

So in the end we get "interconnected" models of the form.
This is important for the definition of our perturbation ports later on.

```asciiart
               ╭────────────────┬───────────╮
               │        voltage │ u out     │
┏━━━━━━━━━━━━━━▽━━┓   ╔═════════△════════╗  │   ╔═══════════════╗
┃ Rest of Network ┃   ║ Shunt-Model      ║  ╰───▷ Device-Model  ║
┃ - voltage in    ┃   ║ - current sum in ║      ║ - voltage in  ║
┃ - current out   ┃   ║ - voltage out    ║  ╭───◁ - current out ║
┗━━━━━━━━━━━━━━▽━━┛   ╚═════════△════════╝  │   ╚═══════════════╝
       current │ i out        ╭─┴─╮         │
               ╰──────────────▷ + ◁─────────╯
                              ╰───╯
```

With that we can start defining the actual models.
For each model we define two buses:
- the shunt bus (dynamic shunt model, static shunt model for powerflow),
- the device bus (machine or inverter model, Slack, PV or PQ model for powerflow) and
- the loopback connection connecting the current-source device to the shunt bus.

### Bus 1: Synchronous Generator (Slack)
=#
ω0 = 2π*50
sg1_bus, bus1, loop1 = let
    @named sm = SyncMachineStatorDynamics(J=3.5, D=1, wL=0.05, R=0.01, ω0=ω0)
    @named sg1_bus = compile_bus(MTKBus(sm); current_source=true)
    set_pfmodel!(sg1_bus, pfSlack(V=1, δ=0; current_source=true, assume_io_coupling=true))

    @named shunt = DynamicRCShunt(R=1/0.6, C=1e-5, ω0=ω0)
    @named bus1 = compile_bus(MTKBus(shunt))
    set_pfmodel!(bus1, pfShunt(G=0.6, B=1e-5))

    loop1 = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:sg1_bus, dst=:bus1)

    sg1_bus, bus1, loop1
end
nothing #hide #md

#=
### Bus 2: Synchronous Generator (PV)
=#
sg2_bus, bus2, loop2 = let
    @named sm = SyncMachineStatorDynamics(J=3.5, D=1, wL=0.05, R=0.01, ω0=ω0)
    @named sg2_bus = compile_bus(MTKBus(sm); current_source=true)
    set_pfmodel!(sg2_bus, pfPV(P=0.5, V=1; current_source=true, assume_io_coupling=true))

    @named shunt = DynamicRCShunt(R=1/0.6, C=1e-5, ω0=ω0)
    @named bus2 = compile_bus(MTKBus(shunt))
    set_pfmodel!(bus2, pfShunt(G=0.6, B=1e-5))

    loop2 = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:sg2_bus, dst=:bus2)

    sg2_bus, bus2, loop2
end
nothing #hide #md

#=
### Bus 3: Grid-Forming Inverter (Droop + LCL)

The GFM uses a [`DroopInverter`](@ref ComposableInverter.DroopInverter) with LCL filter with
cascaded PI controllers for current and voltage in the filter.
The parameters are converted from SimplusGT bandwidth conventions: frequency bandwidths
(e.g. `xfidq=600` Hz) are mapped to PI gains, and cross-coupling feedforward is disabled
(`Fcoupl=0`) to match the MATLAB reference.
=#
gfm_bus, bus3, loop3 = let
    xwLf=0.05; Rf=0.01; xwCf=0.02; xwLc=0.01; Rc=0.002
    Xov=0.01; xDw=0.05; xfdroop=5; xfvdq=300; xfidq=600

    @named droop = ComposableInverter.DroopInverter(;
        filter_type = :LCL,
        droop₊Qset = 0,
        droop₊Kp = xDw*ω0,
        droop₊ω0 = ω0,
        droop₊Kq = 0,
        droop₊τ_q = Inf,
        droop₊τ_p = 1/(xfdroop*2*pi),
        vsrc₊CC1_F = 0,
        vsrc₊CC1_KI = (xfidq*2*pi)^2*(xwLf/ω0)/4,
        vsrc₊CC1_KP = (xfidq*2*pi)*(xwLf/ω0),
        vsrc₊CC1_Fcoupl = 0,
        vsrc₊VC_KP = (xfvdq*2*pi)*(xwCf/ω0),
        vsrc₊VC_KI = (xfvdq*2*pi)^2*(xwCf/ω0)/4*50,
        vsrc₊VC_F = 0,
        vsrc₊VC_Fcoupl = 0,
        vsrc₊X_virt = Xov,
        vsrc₊R_virt = 0,
        vsrc₊Lg = xwLc,
        vsrc₊C = xwCf,
        vsrc₊Rf = Rf,
        vsrc₊Lf = xwLf,
        vsrc₊ω0 = ω0,
        vsrc₊Rg = Rc,
    )
    @named gfm_bus = compile_bus(MTKBus(droop); current_source=true)
    set_pfmodel!(gfm_bus, pfPV(P=0.5, V=1; current_source=true, assume_io_coupling=true))

    @named shunt = DynamicRCShunt(R=1/0.75, C=1e-5, ω0=ω0)
    @named bus3 = compile_bus(MTKBus(shunt))
    set_pfmodel!(bus3, pfShunt(G=0.75, B=1e-5))

    loop3 = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:gfm_bus, dst=:bus3)

    gfm_bus, bus3, loop3
end
nothing #hide #md

#=
### Bus 4: Grid-Following Inverter (DC-link)

The GFL uses a [`SimpleGFLDC`](@ref ComposableInverter.SimpleGFLDC) model with DC-link
dynamics, L filter, PLL with low-pass filter, and current controller.
=#
gfl_bus, bus4, loop4 = let
    V_dc=2.5; C_dc=1.25; f_v_dc=5
    xwLf=0.03; Rf=0.01
    f_pll=5; f_tau_pll=300; f_i_dq=600

    @named gfl = ComposableInverter.SimpleGFLDC(;
        ω0 = ω0,
        Lf = xwLf,
        Rf = Rf,
        PLL_Kp = f_pll*2*pi,
        PLL_Ki = (f_pll*2*pi)^2/4,
        PLL_τ_lpf = 1/(f_tau_pll*2*pi),
        CC1_KP = (xwLf/ω0) * (f_i_dq*2*pi),
        CC1_KI = (xwLf/ω0) * (f_i_dq*2*pi)^2 / 4,
        CC1_F = 0,
        CC1_Fcoupl = 0,
        C_dc = C_dc,
        V_dc = V_dc,
        kp_v_dc = V_dc*C_dc*(f_v_dc*2*pi),
        ki_v_dc = V_dc*C_dc*(f_v_dc*2*pi) * (f_v_dc*2*pi)/4,
    )
    @named gfl_bus = compile_bus(MTKBus(gfl); current_source=true)
    set_pfmodel!(gfl_bus, pfPQ(P=0.5, Q=-0.2; current_source=true))

    @named shunt = DynamicRCShunt(R=1/0.05, C=1e-5, ω0=ω0)
    @named bus4 = compile_bus(MTKBus(shunt))
    set_pfmodel!(bus4, pfShunt(G=0.05, B=1e-5))

    loop4 = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:gfl_bus, dst=:bus4)

    gfl_bus, bus4, loop4
end
nothing #hide #md

#=
## Transmission Lines

All lines use [`DynamicRLBranch`](@ref) (dynamic RL in the rotating dq-frame).
Line 3→4 includes a turns ratio of 0.99. Static [`PiLine`](@ref) models are attached
for the power flow solver.
=#
line12 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=ω0)
    lm = compile_line(MTKLine(branch); name=:l12, src=:bus1, dst=:bus2)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l12_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line23 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=ω0)
    lm = compile_line(MTKLine(branch); name=:l23, src=:bus2, dst=:bus3)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l23_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line31 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=ω0)
    lm = compile_line(MTKLine(branch); name=:l31, src=:bus3, dst=:bus1)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l31_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line34 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=ω0, r_dst=0.99)
    lm = compile_line(MTKLine(branch); name=:l34, src=:bus3, dst=:bus4)
    @named branch_pf = PiLine(R=0.01, X=0.3, r_dst=0.99)
    pfmod = compile_line(MTKLine(branch_pf); name=:l34_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end
nothing #hide #md

#=
## Network Assembly and Initialization

We assemble the network from 8 vertex models (4 device buses + 4 network buses),
4 loopback connections, and 4 transmission lines. Power flow is solved and used to
initialize all dynamic states.

After defining the network we solve the powerflow and initialize the dynamic states
at the powerflow solution. This is the operating point around which we will analyze the system.
=#
nw = Network([sg1_bus, bus1, sg2_bus, bus2, gfm_bus, bus3, gfl_bus, bus4],
    [loop1, loop2, loop3, loop4, line12, line23, line31, line34]; warn_order=false)

pfs = solve_powerflow(nw; abstol=1e-10, reltol=1e-10)
s0 = initialize_from_pf!(nw; pfs, tol=1e-7, nwtol=1e-7)
nothing #hide #md

#=
## Eigenvalue Analysis

[`jacobian_eigenvals`](@extref NetworkDynamics.jacobian_eigenvals) linearizes the system,
eliminates algebraic constraints via Schur complement, and returns the eigenvalues of the
reduced state matrix. We divide by ``2\pi`` to convert from rad/s to Hz.
=#
eigenvalues = jacobian_eigenvals(nw, s0) ./ (2 * pi)

let
    fig = Figure(size=(600,500))
    ax1 = Axis(fig[1, 1], xlabel="Real Part [Hz]", ylabel="Imaginary Part [Hz]", title="Global Pole Map")
    scatter!(ax1, real.(eigenvalues), imag.(eigenvalues), marker=:xcross)
    ax2 = Axis(fig[1, 2], xlabel="Real Part [Hz]", ylabel="Imaginary Part [Hz]", title="Zoomed In")
    scatter!(ax2, real.(eigenvalues), imag.(eigenvalues), marker=:xcross)
    xlims!(ax2, -80, 20); ylims!(ax2, -150, 150)
    fig
end

#=
The results match what we get from SimplusGT.

!!! details "SimplusGT Reference: Pole Map"
    ![image](../assets/SimplusGTPlots/Figure_100.png)
=#

#=
## Impedance-Based Bode Analysis

SimplusGT models each device as a transfer function from bus voltage to injected current:
the device admittance $G(s)$. To get a closed loop model, they add feedback through the
network impedance matrix $Z(s)$.
For impedance-based stability analysis, they directionally perturb the voltage input of a
device *without* perturbing that voltage point in the network -- effectively opening the
loop at that bus.

```asciiart
             (Admittance-like)
            ╭────────────────╮
δu_dq╶─→●──→┤ Device TF G(s) ├→───●──╴ i_dq
        │   ╰────────────────╯    │
        │                         │
        │   ╭────────────────╮    │
(bus    ╰──←┤  Grid TF Z(s)  ├←───╯
 voltage)   ╰────────────────╯
             (Impedance-like)
```

This kind of analysis can be replicated using the [`linearize_network`](@extref NetworkDynamics.linearize_network-Tuple{})
function while proviing the `in` and `out` keyword aguments.
With `in` it is possible to specify the perturbation port, `out` defines the observed stated under this output.

Perturbation in the linearization sense is *directed*, our device models however are fully acausal, there is
no directionality within their definition.
Therefore we a re bit restricted by the choice of perturbation channel: it has to be a "connection" between network
components, i.e. either a vertex input, a vetex output, an edge input or an edge output.

Considering the component interconnection on the network level sketched out above, we have 6 channels
for directional perturbation injection, 4 of which are unique:
```asciiart
               ╭────────────────┬───────────╮
        δu_in →∙        δu_out →∙           ∙← δu_in
        (line) │         (hub)  │           │ (device)
┏━━━━━━━━━━━━━━▽━━┓   ╔═════════△════════╗  │   ╔═══════════════╗
┃ Rest of Network ┃   ║ Shunt-Model      ║  ╰───▷ Device-Model  ║
┃ - voltage in    ┃   ║ - current sum in ║      ║ - voltage in  ║
┃ - current out   ┃   ║ - voltage out    ║  ╭───◁ - current out ║
┗━━━━━━━━━━━━━━▽━━┛   ╚═════════△════════╝  │   ╚═══════════════╝
       δi_out →∙        δi_out →∙           ∙← δi_out
       (line)  │        (hub) ╭─┴─╮         │ (device)
               ╰──────────────▷ + ◁─────────╯
                              ╰───╯
```
- `δu_in (line)`: Perturbation of the voltage at the line side (perturbs line)
- `δu_out (hub)`: Perturbation of the voltage at the hub output (influences both device and line)
- `δu_in (device)`: Perturbation of the voltage at the device input (perturbs device)
- `δi_out (line)`, `δi_out (hub)` and `δi_out (device)`: Perturbation of the current before or after aggregation. Mathematicially equivalent, only the hub input either way.

First, we are interested in the $Y_{dd}(s)$ admittance, which is the transfer function from `δu_in (device)` to the device output.
I.e. we perturb the input voltage for the device model and observe the change in current output. This is directly
equivalent to the ports used in SimplusGT shown above.

Since we are only interested in $Y_dd$ for now, we get away with a SISO system. We just specify a single perturbance
and single obsered channel:
=#
(; M, A, B, C, D, G) = linearize_network(nw, s0; in=VIndex(:sg1_bus, :busbar₊u_r), out=VIndex(:sg1_bus, :busbar₊i_r))
nothing #hide
#=
We get back a named tuple representing our descriptor system. It contains all matrices as well as a `G(s)` transfer function.
You can use the matrices to construct systems in ControlSystems.jl or use other libraries to work with them.
Here, we are only intersted in the transfer function. Since we can evaluate that directly:
=#
G(im*2π*50)
#=
We can define our own bode plotting function without depending on additional external libraries:
=#
function bode_plot(Gs, title="", labels=["Bus $i" for i in 1:length(Gs)])
    ## make colors match the matlab plot better
    with_theme(Theme(palette = (; color = Makie.wong_colors()[[1, 6, 2, 4]]))) do
        fig = Figure(; size=(800, 600))
        Label(fig[1, 1], title*"Bode Plot", fontsize=16, halign=:center, tellwidth=false)
        ax1 = Axis(fig[2, 1], xlabel="Frequency (rad/s)", ylabel="Gain (dB)", xscale=log10)
        ax2 = Axis(fig[3, 1], xlabel="Frequency (rad/s)", ylabel="Phase (deg)", xscale=log10)
        for (G, label) in zip(Gs, labels)
            fs = 10 .^ (range(log10(1e-1), log10(1e4); length=1000))
            ωs = 2π * fs
            ss = im .* ωs
            gains = map(s -> 20 * log10(abs(G(s))), ss)
            phases = rad2deg.(unwrap_rad(map(s -> angle(G(s)), ss)))
            lines!(ax1, fs, gains; label, linewidth=2)
            lines!(ax2, fs, phases; label, linewidth=2)
        end
        axislegend(ax1)
        fig
    end
end
nothing #hide #md
#=
In order to reproduce the $Y_{dd}$ bode plot, we just need to generate the linearization around
every component and plot the results.
=#
Gs = map([:sg1_bus, :sg2_bus, :gfm_bus, :gfl_bus]) do COMP
    vs = VIndex(COMP, :busbar₊u_r)
    cs = VIndex(COMP, :busbar₊i_r)
    G = NetworkDynamics.linearize_network(nw, s0; in=vs, out=cs).G
end
bode_plot(Gs, "Y_dd ")

#=
There is a small difference in unwrapping algorithm, but besides that we nicely replicate the results from SimplusGT.

!!! details "SimplusGT Reference: Bode Plot"
    ![image](../assets/SimplusGTPlots/Figure_200.png)
=#

#=
## Time-Domain Simulation

Since we have the full nonlinear model at hand, lets also do a time domain simualtion.
To perturb the system, we apply a three-phase short circuit at Bus 1 by reducing the shunt resistance to near zero for a short time.
The faul starts a 0.1s and is cleared at 0.2s.
=#
affect = ComponentAffect([], [:shunt₊R]) do u, p, ctx
    if ctx.t == 0.1
        println("Short Circuit at Bus 1 at t=0.1s")
        p[:shunt₊R] = 1e-6
    elseif ctx.t == 0.2
        println("Clearing Short Circuit at Bus 1 at t=0.2s")
        p[:shunt₊R] = 1/0.6
    end
end
short = PresetTimeComponentCallback([0.1, 0.2], affect)
prob = ODEProblem(nw, s0, (0,30); add_comp_cb=VIndex(:bus1)=>short)
sol = solve(prob, Rodas5P())
nothing # hide

#=
Lets compare the results against a similar EMT simulation done in the Simulink model generated by SimplusGT.
=#
with_theme(Theme(palette = (; color = Makie.wong_colors()[[1, 6, 2, 4]]))) do
    fig = Figure()
    ax = Axis(fig[1,1], title="Voltage Magnitude", limits=((0,0.5), (0.0, 1.2)))
    ts = refine_timeseries(sol.t)
    lines!(ax, ts, sol(ts, idxs=VIndex(:bus1, :busbar₊u_mag)).u; label="Bus 1")
    lines!(ax, ts, sol(ts, idxs=VIndex(:bus2, :busbar₊u_mag)).u; label="Bus 2")
    lines!(ax, ts, sol(ts, idxs=VIndex(:bus3, :busbar₊u_mag)).u; label="Bus 3")
    lines!(ax, ts, sol(ts, idxs=VIndex(:bus4, :busbar₊u_mag)).u; label="Bus 4")
    ax2 = Axis(fig[2,1], limits=((0,30), (0.0, 1.2)))
    lines!(ax2, ts, sol(ts, idxs=VIndex(:bus1, :busbar₊u_mag)).u; label="Bus 1")
    lines!(ax2, ts, sol(ts, idxs=VIndex(:bus2, :busbar₊u_mag)).u; label="Bus 2")
    lines!(ax2, ts, sol(ts, idxs=VIndex(:bus3, :busbar₊u_mag)).u; label="Bus 3")
    lines!(ax2, ts, sol(ts, idxs=VIndex(:bus4, :busbar₊u_mag)).u; label="Bus 4")
    axislegend(ax2; position=:rb)
    fig
end

#=
!!! details "SimplusGT Reference: Short Circuit Response"
    ![image](../assets/SimplusGTPlots/Voltage_Trajectory.png)
=#


