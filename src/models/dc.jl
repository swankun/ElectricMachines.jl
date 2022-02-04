export DCMotor, simulate

abstract type AbstractDCMotor end

function defaultparams()
    return (
        R = 2,
        L = 0.002,
        Kb = 0.07,
        J = 6e-5,
        f = 0.4e-3,
        KT = 0.07,
        ppr = 500,
        curret_gains = (kp=1e2, ki=1e4),
        observer_poles = (p1=50., p2=200., p3=10.),
        Vmax = 48,
        imax = 10,
    )
end

struct DCMotor{F,PS,IREF,TAU} <: AbstractDCMotor
    eom!::F
    params::PS
    iref::IREF
    load::TAU
end

# function DCMotor(; iref::Function=zero, load::Function=zero, init_params::Function=defaultparams)
#     params = init_params()
#     R, L, Kb, J, f, KT, _, current_gains, observer_poles, Vmax, _ = params
#     p1, p2, p3 = observer_poles
#     Kp, Ki = current_gains
#     ℓ1 = p1 + p2 + p3 - f/J
#     ℓ2 = p1*p2 + p1*p3 + p2*p3 - ℓ1*f/J
#     ℓ3 = -p1*p2*p3
#     f!(dx,x,p,t) = begin
#         i, ω, θ, i_error, θbar, ωbar, τbar, encoder = x
#         p_error = iref(t) - i
#         effort = Kp*p_error + Ki*i_error
#         v = clamp(effort, -Vmax, Vmax)
#         dx[1] = (-R*i - Kb*ωbar + v) / L
#         dx[2] = (-f*ωbar + KT*i - load(t)) / J
#         dx[3] = ω
#         dx[4] = p_error
#         dx[5] = ωbar + ℓ1*(encoder-θbar)
#         dx[6] = (-f*ωbar + KT*i - τbar) / J + ℓ2*(encoder-θbar)
#         dx[7] = ℓ3*(encoder-θbar)*J
#         dx[8] = zero(encoder)
#         nothing
#     end
#     DCMotor{typeof(f!),typeof(params),typeof(iref),typeof(load)}(f!, params, iref, load)
# end

function DCMotor(; iref::Function=zero, load::Function=zero, init_params::Function=defaultparams)
    θref, ωref = homework4()
    params = init_params()
    R, L, Kb, J, f, KT, _, current_gains, observer_poles, Vmax, _ = params
    Kp, Ki = current_gains
    r1, r2, r3 = (10., 20., 30.)
    k0 = r1*r2*r3
    k1 = r1*r2 + r1*r3 + r2*r3
    k2 = r1 + r2 + r3 - f/J
    p1, p2, p3 = observer_poles
    ℓ1 = p1 + p2 + p3 - f/J
    ℓ2 = p1*p2 + p1*p3 + p2*p3 - ℓ1*f/J
    ℓ3 = -p1*p2*p3
    f!(dx,x,p,t) = begin
        i, ω, θ, i_ierror, θbar, ωbar, τbar, encoder, pos_ierror = x
        pos_perror = θref(t) - encoder
        vel_perror = ωref(t) - ωbar
        w =  -(k0*pos_ierror + k1*pos_perror + k2*vel_perror)
        ir = iref(t) - J/KT*w
        i_perror = ir - i
        effort = Kp*i_perror + Ki*i_ierror
        # v = clamp(effort, -Vmax, Vmax)
        v = effort
        dx[1] = (-R*i - Kb*ωbar + v) / L
        dx[2] = (-f*ωbar + KT*i - load(t)) / J
        dx[3] = ω
        dx[4] = i_perror
        dx[5] = ωbar + ℓ1*(encoder-θbar)
        dx[6] = (-f*ωbar + KT*i - τbar) / J + ℓ2*(encoder-θbar)
        dx[7] = ℓ3*(encoder-θbar)*J
        dx[8] = zero(encoder)
        dx[9] = pos_perror
        nothing
    end
    DCMotor{typeof(f!),typeof(params),typeof(iref),typeof(load)}(f!, params, iref, load)
end

function simulate(dc::DCMotor; Δt=1e-3, tf=1.0, kwargs...)
    prob = ODEProblem{true}(dc.eom!, zeros(9), (zero(tf),tf))
    N = 2π/dc.params.ppr/4
    affect!(integrator) = begin
        integrator.u[8] = div(integrator.u[3], N)*N
    end
    cb = PeriodicCallback(affect!, Δt)
    return solve(prob, Tsit5(), callback=cb; kwargs...)
end

