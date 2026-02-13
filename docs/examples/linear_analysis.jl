#=
# [Linear Analysis of a 4-Bus System](@id linear-analysis)

This example can be downloaded as a normal Julia script [here](@__NAME__.jl). #md

This example demonstrates **eigenvalue analysis** and **impedance-based Bode analysis** of a
mixed-technology 4-bus power system. The system parameters and topology are taken from
[SimplusGT](https://github.com/Future-Power-Networks/Simplus-Grid-Tool) (Simplus Grid Tool),
and eigenvalue results have been validated against the MATLAB reference.

Please note that this tutorial only replicates a small part of what SimplusGT is capable of,
so please check out this greate toolbox for in-depth linear analysis of powersystems!

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
nothing #hide

#=
## Synchronous Machine with Stator Dynamics
In order to replicate the results we need to use identical or nearly identical models to
SimplusGT. For that we need to implement a different kind of synchronous machine that
- retains the stator flux dynamics (i.e. defines the current output as an differential equation),
- assumes constant field flux (parameter which is initialized at operation point) and
- shows no governor dynamics (i.e. mechanical torque is a constant parameter initialized at operation point).

The model matches the equations of the Type 0 Apparatus in SimplusGT.
=#
@mtkmodel SyncMachineStatorDynamics begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        J, [description="Inertia constant [MWs²/MVA]"]
        D, [description="Damping coefficient [pu]"]
        wL, [description="Stator inductance * base frequency [pu]"]
        R, [description="Stator resistance [pu]"]
        w0, [description="Base frequency [rad/s]"]
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
        J_pu = J*2/w0^2
        D_pu = D/w0^2
        L = wL/w0
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
        Dt(theta) ~ w - w0
    end
end
nothing #hide #md

#=
## 4-Bus Network Setup

The network replicates the default 4-bus example from SimplusGT's `UserData.json`:

| Bus | Device | PF Type | Key Parameters |
|:----|:-------|:--------|:---------------|
| 1 | SG Type 0 | Slack (V=1) | J=3.5, D=1, wL=0.05, R=0.01 |
| 2 | SG Type 0 | PV (P=0.5, V=1) | J=3.5, D=1, wL=0.05, R=0.01 |
| 3 | GFM Type 20 (Droop, LCL) | PV (P=0.5, V=1) | Droop + voltage/current control |
| 4 | GFL Type 10 (DC-link) | PQ (P=0.5, Q=-0.2) | PLL + DC voltage control |

Lines 1↔2, 2↔3, 3↔1 have R=0.01, wL=0.3. Line 3→4 adds a turns ratio of 0.99.

!!! note EMT Models
    The entire modeling in SimplusGT uses **EMT Components**, meaning that there is no
    static dq-phasor calculation like `i_dq = Z_dq⋅u_dq`. Instead, every   mod

### Current Source Buses and Loopback Connections

SimplusGT models each device as a transfer function from bus voltage to injected current
(an admittance ``Y(s)``), and closes the loop through the network admittance matrix.

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



```
╭─────────╮                       ╭──╴Powerline 1
│ Current │          Vbus         │
│ Source  ├────→──────●────←──────●
│ Device  │ i_device  ┴  i_grid   │
╰─────────╯           ┬           ╰──╴Powerline 2
                      ⏚

```

For
impedance-based stability analysis, they directionally perturb the voltage input of a
device *without* simultaneously perturbing that voltage point in the network -- effectively
opening the loop at that bus.

We replicate this using the current-source bus + loopback pattern: each device bus
(`current_source=true`) is paired with a network bus carrying a `DynamicRCShunt` via a
[`LoopbackConnection`](@extref NetworkDynamics.LoopbackConnection). This gives us separate
voltage/current access points that can be independently targeted by
[`linearize_network`](@extref NetworkDynamics.linearize_network). See [Injector
Nodes](@extref injector-nodes) in the NetworkDynamics.jl docs for background on this pattern.
=#

w0 = 2π*50

#=
### Bus 1: Synchronous Generator (Slack)
=#
sg1_bus, bus1, loop1 = let
    @named sm = SyncMachineStatorDynamics(J=3.5, D=1, wL=0.05, R=0.01, w0=w0)
    @named sg1_bus = compile_bus(MTKBus(sm); current_source=true)
    set_pfmodel!(sg1_bus, pfSlack(V=1, δ=0; current_source=true, assume_io_coupling=true))

    @named shunt = DynamicRCShunt(R=1/0.6, C=1e-5, ω0=w0)
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
    @named sm = SyncMachineStatorDynamics(J=3.5, D=1, wL=0.05, R=0.01, w0=w0)
    @named sg2_bus = compile_bus(MTKBus(sm); current_source=true)
    set_pfmodel!(sg2_bus, pfPV(P=0.5, V=1; current_source=true, assume_io_coupling=true))

    @named shunt = DynamicRCShunt(R=1/0.6, C=1e-5, ω0=w0)
    @named bus2 = compile_bus(MTKBus(shunt))
    set_pfmodel!(bus2, pfShunt(G=0.6, B=1e-5))

    loop2 = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:sg2_bus, dst=:bus2)

    sg2_bus, bus2, loop2
end
nothing #hide #md

#=
### Bus 3: Grid-Forming Inverter (Droop + LCL)

The GFM uses a [`DroopInverter`](@ref ComposableInverter.DroopInverter) with LCL filter.
Parameters are converted from SimplusGT bandwidth conventions: frequency bandwidths
(e.g. `xfidq=600` Hz) are mapped to PI gains, and cross-coupling feedforward is disabled
(`Fcoupl=0`) to match the MATLAB reference.
=#
gfm_bus, bus3, loop3 = let
    xwLf=0.05; Rf=0.01; xwCf=0.02; xwLc=0.01; Rc=0.002
    Xov=0.01; xDw=0.05; xfdroop=5; xfvdq=300; xfidq=600

    @named droop = ComposableInverter.DroopInverter(;
        filter_type=:LCL,
        droop₊Qset = 0,
        droop₊Kp = xDw*w0,
        droop₊ω0 = w0,
        droop₊Kq = 0,
        droop₊τ_q = Inf,
        droop₊τ_p = 1/(xfdroop*2*pi),
        vsrc₊CC1_F = 0,
        vsrc₊CC1_KI = (xfidq*2*pi)^2*(xwLf/w0)/4,
        vsrc₊CC1_KP = (xfidq*2*pi)*(xwLf/w0),
        vsrc₊CC1_Fcoupl = 0,
        vsrc₊VC_KP = (xfvdq*2*pi)*(xwCf/w0),
        vsrc₊VC_KI = (xfvdq*2*pi)^2*(xwCf/w0)/4*50,
        vsrc₊VC_F = 0,
        vsrc₊VC_Fcoupl = 0,
        vsrc₊X_virt = Xov,
        vsrc₊R_virt = 0,
        vsrc₊Lg = xwLc,
        vsrc₊C = xwCf,
        vsrc₊Rf = Rf,
        vsrc₊Lf = xwLf,
        vsrc₊ω0 = w0,
        vsrc₊Rg = Rc,
    )
    @named gfm_bus = compile_bus(MTKBus(droop); current_source=true)
    set_pfmodel!(gfm_bus, pfPV(P=0.5, V=1; current_source=true, assume_io_coupling=true))

    @named shunt = DynamicRCShunt(R=1/0.75, C=1e-5, ω0=w0)
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
        ω0 = w0,
        Lf = xwLf,
        Rf = Rf,
        PLL_Kp = f_pll*2*pi,
        PLL_Ki = (f_pll*2*pi)^2/4,
        PLL_τ_lpf = 1/(f_tau_pll*2*pi),
        CC1_KP = (xwLf/w0) * (f_i_dq*2*pi),
        CC1_KI = (xwLf/w0) * (f_i_dq*2*pi)^2 / 4,
        CC1_F = 0,
        CC1_Fcoupl = 0,
        C_dc = C_dc,
        V_dc = V_dc,
        kp_v_dc = V_dc*C_dc*(f_v_dc*2*pi),
        ki_v_dc = V_dc*C_dc*(f_v_dc*2*pi) * (f_v_dc*2*pi)/4,
    )
    @named gfl_bus = compile_bus(MTKBus(gfl); current_source=true)
    set_pfmodel!(gfl_bus, pfPQ(P=0.5, Q=-0.2; current_source=true))

    @named shunt = DynamicRCShunt(R=1/0.05, C=1e-5, ω0=w0)
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
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=w0)
    lm = compile_line(MTKLine(branch); name=:l12, src=:bus1, dst=:bus2)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l12_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line23 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=w0)
    lm = compile_line(MTKLine(branch); name=:l23, src=:bus2, dst=:bus3)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l23_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line31 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=w0)
    lm = compile_line(MTKLine(branch); name=:l31, src=:bus3, dst=:bus1)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l31_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line34 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=w0, r_dst=0.99)
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
=#
nw = Network([sg1_bus, bus1, sg2_bus, bus2, gfm_bus, bus3, gfl_bus, bus4],
    [loop1, loop2, loop3, loop4, line12, line23, line31, line34]; warn_order=false)

pfs = solve_powerflow(nw; abstol=1e-10, reltol=1e-10)
show_powerflow(pfs)
s0 = initialize_from_pf!(nw; pfs, subverbose=true, tol=1e-7, nwtol=1e-7)
nothing #hide #md

#=
## Eigenvalue Analysis

[`jacobian_eigenvals`](@extref NetworkDynamics.jacobian_eigenvals) linearizes the system,
eliminates algebraic constraints via Schur complement, and returns the eigenvalues of the
reduced state matrix. We divide by ``2\pi`` to convert from rad/s to Hz.
=#
eigenvalues = jacobian_eigenvals(nw, s0) ./ (2 * pi)
println("4-Bus Eigenvalues:")
display(eigenvalues)

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
!!! details "SimplusGT Reference: Pole Map"
    *Reference plot not yet exported. Will be added in a future update.*
=#

#=
## Impedance-Based Bode Analysis

SimplusGT models each device as a transfer function from bus voltage to injected current
(an admittance ``Y(s)``), and closes the loop through the network admittance matrix.
For impedance-based stability analysis, they directionally perturb the voltage input of a
device *without* perturbing that voltage point in the network -- effectively opening the
loop at that bus.

We replicate this using `linearize_network` with `in` and `out` keyword arguments. By
specifying the voltage states of a device bus as inputs and the current states as outputs,
we obtain the device admittance transfer function ``Y_{dd}(s)``. The loopback connection
topology naturally provides the separation: perturbing the voltage input of the device bus
does not simultaneously perturb the network bus voltage.
=#

function bode_plot(Gs, title="", labels=["Bus $i" for i in 1:length(Gs)])
    with_theme(Theme(palette = (; color = Makie.wong_colors()[[1, 6, 2, 4]]))) do
        fig = Figure(; size=(800, 600))
        Label(fig[1, 1], title*"Bode Plot", fontsize=16, halign=:center, tellwidth=false)
        ax1 = Axis(fig[2, 1], xlabel="Frequency (rad/s)", ylabel="Gain (dB)", xscale=log10)
        ax2 = Axis(fig[3, 1], xlabel="Frequency (rad/s)", ylabel="Phase (deg)", xscale=log10)

        for (G, label) in zip(Gs, labels)
            fs = 10 .^ (range(log10(1e-1), log10(1e4); length=1000))
            ωs = 2π * fs
            ss = im .* ωs
            output, input = 1, 1
            gains = map(s -> 20 * log10(abs(G(s)[output, input])), ss)
            phases = rad2deg.(unwrap_rad(map(s -> angle(G(s)[output, input]), ss)))

            lines!(ax1, fs, gains; label, linewidth=2)
            lines!(ax2, fs, phases; label, linewidth=2)
        end
        axislegend(ax1)
        fig
    end
end
nothing #hide #md

# Each device's admittance is obtained by linearizing with voltage as input and current as output.
Gs = map([:sg1_bus, :sg2_bus, :gfm_bus, :gfl_bus]) do COMP
    vs = VIndex(COMP, [:busbar₊u_r, :busbar₊u_i])
    cs = VIndex(COMP, [:busbar₊i_r, :busbar₊i_i])
    G = NetworkDynamics.linearize_network(nw, s0; in=vs, out=cs).G
end
bode_plot(Gs, "Y_dd ")

#=
!!! details "SimplusGT Reference: Bode Plot"
    *Reference plot not yet exported. Will be added in a future update.*
=#

#=
## Time-Domain Simulation

We validate the linearization results against a nonlinear time-domain simulation.
A three-phase short circuit is applied at Bus 1 (shunt resistance reduced to near-zero)
at t=0.1s and cleared at t=0.2s.
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

# ### Voltage Response
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
    *Reference plot not yet exported. Will be added in a future update.*
=#
