export homework4

function homework4()
    T = 70e-2
    xf = π
    t1 = 3T/10
    t2 = 7T/10
    ωmax = xf/t2
    c = [t1^2 t1^3; 2t1 3t1^2] \ [ωmax; 0]
    θ(t) = begin
        if 0 ≤ t ≤ t1
            return c[1]/3*t^3 + c[2]/4*t^4
        elseif t1 ≤ t ≤ t2
            return ωmax*t1/2 + ωmax*(t-t1)
        elseif t2 ≤ t ≤ T
            return ωmax*t2 - c[1]/3*(T-t)^3 - c[2]/4*(T-t)^4
        else
            return xf
        end
    end
    ω(t) = begin
        if 0 ≤ t ≤ t1
            return c[1]*t^2 + c[2]*t^3
        elseif t1 ≤ t ≤ t2
            return ωmax
        elseif t2 ≤ t ≤ T
            return c[1]*(T-t)^2 + c[2]*(T-t)^3
        else
            return zero(t)
        end
    end
    α(t) = begin
        if 0 ≤ t ≤ t1
            return 2c[1]*t + 3c[2]*t^2
        elseif t1 ≤ t ≤ t2
            return zero(t)
        elseif t2 ≤ t ≤ T
            return -2c[1]*(T-t) - 3c[2]*(T-t)^2
        else
            return zero(t)
        end
    end
    p = defaultparams()
    i(t) = (p.J*α(t) + p.f*ω(t))/p.KT
    didt(t) = ForwardDiff.derivative(i,t)
    v(t) = p.L*didt(t) + p.R*i(t) + p.Kb*ω(t)
    fig = Figure()
    t = range(0, T, length=1001)
    # return map(ω, t)
    lines(fig[1,1], t, θ)
    lines(fig[1,2], t, ω)
    lines(fig[2,1], t, i)
    lines(fig[2,2], t, v)
    save("plots/traj.png", fig)
    return θ, ω
end
