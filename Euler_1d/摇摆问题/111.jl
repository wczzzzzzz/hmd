using GLMakie
using LinearAlgebra
using TimerOutputs

g = 9.81
L1, L2 = 1.0, 1.0
m1, m2 = 1.0, 1.0

θ1₀, θ2₀ = 0.5, 0.3
θ̇1₀, θ̇2₀ = 0.0, 0.0

T = 20.0
Δt = 0.01
t_steps = 0:Δt:T
n_t = length(t_steps)

function compute_accelerations(θ1, θ2, θ̇1, θ̇2)
    Δθ = θ1 - θ2
    num_θ1 = -m2*L1*L2*cos(Δθ)*θ̇2^2 * sin(Δθ) - (m1+m2)*g*L1*sin(θ1)
    den_θ1 = (m1+m2)*L1^2 + m2*L1^2*cos(Δθ)^2
    θ̈1 = num_θ1 / den_θ1  

    num_θ2 = m2*L1*L2*cos(Δθ)*θ̇1^2 * sin(Δθ) - m2*g*L2*sin(θ2)
    den_θ2 = m2*L2^2 + m2*L1^2*cos(Δθ)^2
    θ̈2 = num_θ2 / den_θ2  

    return θ̈1, θ̈2
end

function explicit_euler()
    θ1 = zeros(n_t)
    θ2 = zeros(n_t)
    θ̇1 = zeros(n_t)
    θ̇2 = zeros(n_t)
    θ1[1], θ2[1] = θ1₀, θ2₀
    θ̇1[1], θ̇2[1] = θ̇1₀, θ̇2₀

    for n in 1:n_t-1
        θ̈1, θ̈2 = compute_accelerations(θ1[n], θ2[n], θ̇1[n], θ̇2[n])
        
        θ̇1[n+1] = θ̇1[n] + Δt * θ̈1
        θ̇2[n+1] = θ̇2[n] + Δt * θ̈2
        θ1[n+1] = θ1[n] + Δt * θ̇1[n]
        θ2[n+1] = θ2[n] + Δt * θ̇2[n]
    end
    return θ1, θ2, θ̇1, θ̇2
end

function residuals(x, θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt)
    θ1_next, θ2_next, θ̇1_next, θ̇2_next = x
    Δθ = θ1_next - θ2_next

    num_θ1 = -m2*L1*L2*cos(Δθ)*θ̇2_next^2 * sin(Δθ) - (m1+m2)*g*L1*sin(θ1_next)
    den_θ1 = (m1+m2)*L1^2 + m2*L1^2*cos(Δθ)^2 
    θ̈1_next = num_θ1 / den_θ1  

    num_θ2 = m2*L1*L2*cos(Δθ)*θ̇1_next^2 * sin(Δθ) - m2*g*L2*sin(θ2_next)
    den_θ2 = m2*L2^2 + m2*L1^2*cos(Δθ)^2 
    θ̈2_next = num_θ2 / den_θ2  

    r1 = θ1_next - θ1_curr - Δt * θ̇1_next
    r2 = θ2_next - θ2_curr - Δt * θ̇2_next
    r3 = θ̇1_next - θ̇1_curr - Δt * θ̈1_next
    r4 = θ̇2_next - θ̇2_curr - Δt * θ̈2_next
    return [r1, r2, r3, r4]
end

function jacobian(x, θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt, ε=1e-8)
    J = zeros(4, 4)
    f0 = residuals(x, θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt)
    for i in 1:4
        x_eps = copy(x)
        x_eps[i] += ε
        f_eps = residuals(x_eps, θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt)
        J[:, i] = (f_eps - f0) / ε
    end
    return J
end

function implicit_euler()
    θ1 = zeros(n_t)
    θ2 = zeros(n_t)
    θ̇1 = zeros(n_t)
    θ̇2 = zeros(n_t)
    θ1[1], θ2[1] = θ1₀, θ2₀
    θ̇1[1], θ̇2[1] = θ̇1₀, θ̇2₀

    for n in 1:n_t-1
        θ1_curr, θ2_curr = θ1[n], θ2[n]
        θ̇1_curr, θ̇2_curr = θ̇1[n], θ̇2[n]

        θ̈1_guess, θ̈2_guess = compute_accelerations(θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr)
        θ̇1_guess = θ̇1_curr + Δt * θ̈1_guess
        θ̇2_guess = θ̇2_curr + Δt * θ̈2_guess
        θ1_guess = θ1_curr + Δt * θ̇1_guess
        θ2_guess = θ2_curr + Δt * θ̇2_guess
        x = [θ1_guess, θ2_guess, θ̇1_guess, θ̇2_guess]

        for iter in 1:20
            r = residuals(x, θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt)
            if norm(r) < 1e-8
                break
            end
            J = jacobian(x, θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt)
            x -= J \ r
        end

        θ1[n+1], θ2[n+1], θ̇1[n+1], θ̇2[n+1] = x
    end
    return θ1, θ2, θ̇1, θ̇2
end

to = TimerOutput()
@timeit to "显式欧拉" θ1_exp, θ2_exp, θ̇1_exp, θ̇2_exp = explicit_euler()
@timeit to "隐式欧拉" θ1_imp, θ2_imp, θ̇1_imp, θ̇2_imp = implicit_euler()


fig = Figure(resolution=(1200, 600))
ax_exp = Axis(fig[1,1], title="显式欧拉", xlabel="时间 (s)", ylabel="角度 (rad)")
ax_imp = Axis(fig[1,2], title="隐式欧拉", xlabel="时间 (s)", ylabel="角度 (rad)")

lines!(ax_exp, t_steps, θ1_exp, color=:blue, label="θ₁")
lines!(ax_exp, t_steps, θ2_exp, color=:red, label="θ₂")
lines!(ax_imp, t_steps, θ1_imp, color=:blue, label="θ₁")
lines!(ax_imp, t_steps, θ2_imp, color=:red, label="θ₂")

axislegend(ax_exp; position=:rb)
axislegend(ax_imp; position=:rb)
display(fig)

show(to)