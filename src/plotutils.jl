export default_plot_theme, plotobserver, plotfeedback

function default_plot_theme()
    majorfontsize = 36*1.5
    minorfontsize = 24*1.5
    T = Theme(
        Axis = (
            xlabelfont="CMU Serif",
            ylabelfont="CMU Serif",
            xticklabelfont="CMU Serif",
            yticklabelfont="CMU Serif",
            titlefont="CMU Serif",
            xlabelsize=minorfontsize,
            ylabelsize=minorfontsize,
            xticklabelsize=minorfontsize,
            yticklabelsize=minorfontsize,
            titlesize=majorfontsize,
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false,
        ),
        Lines = (
            linewidth = 3,
        ),
        Legend = (
            labelfont = "CMU Serif",
            labelsize = minorfontsize
        )
    )
    set_theme!(T)
end

function plotobserver(
    sys::DCMotor=DCMotor(iref=(t)->0.5sin(10π*t)), 
    fig=Figure(resolution=(1800,600)); 
    kwargs...
)
    sol = simulate(sys; kwargs...)
    t = range(first(sol.t), stop=last(sol.t), length=1001)
    x = Array(sol(t))
    ω = x[2,:]
    ωbar = x[6,:]
    τbar = x[7,:]
    
    ax1 = Axis(fig[1,1],
        xlabel="Time (s)", title="Angular velocity (rad/s)"
    )
    lines!(ax1, t, ω, color=:blue, label=L"\omega(t)")
    lines!(ax1, t, ωbar, color=:red, linedwidth=3, label=L"\hat{\omega}(t)")
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1,2],
        xlabel="Time (s)", title="Load torque (Nm)"
    )
    lines!(ax2, t, sys.load, color=:blue, label=L"\tau_{L}(t)")
    lines!(ax2, t, τbar, color=:red, linedwidth=3, label=L"\hat{\tau}_{L}(t)")
    axislegend(ax2, position=:rb)
    save("plots/out.png", fig)
    nothing
end

function plotfeedback(
    sys::DCMotor=DCMotor(), 
    fig=Figure(resolution=(1800,600)); 
    kwargs...
)
    sol = simulate(sys; kwargs...)
    t = range(first(sol.t), stop=last(sol.t), length=1001)
    x = Array(sol(t))
    ω = x[2,:]
    θ = x[3,:]
    
    ax1 = Axis(fig[1,1],
        xlabel="Time (s)", title="Angular velocity (rad/s)"
    )
    lines!(ax1, t, ω, color=:blue, label=L"\omega(t)")

    ax2 = Axis(fig[1,2],
        xlabel="Time (s)", title="Position (rad)"
    )
    lines!(ax2, t, θ, color=:blue, label=L"\theta(t)")
    save("plots/out.png", fig)
    nothing
end
