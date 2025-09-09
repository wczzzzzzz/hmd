using GLMakie
using LinearAlgebra

g = 9.81 
l1 = 1.0
l2 = 1.0
m1 = 1.0
m2 = 1.0

θ1₀ = 0.5
θ2₀ = 0.3
θ̇1₀ = 0.0
θ̇2₀ = 0.0

T = 20.0  
Δt = 0.01 
t = 0.0:Δt:T
nₜ = length(t)

θ1 = zeros(nₜ)
θ2 = zeros(nₜ)
θ̇1 = zeros(nₜ)
θ̇2 = zeros(nₜ)

θ1[1] = θ1₀
θ2[1] = θ2₀
θ̇1[1] = θ̇1₀
θ̇2[1] = θ̇2₀

function compute_accelerations(θ1, θ2, θ̇1, θ̇2)
    Δθ = θ1 - θ2
    d1 = (m1 + m2)*l1 - m2*l1*cos(Δθ)^2
    d2 = (l2/l1)*d1
    
    θ̈1 = ( -g*(2m1 + m2)*sin(θ1) - m2*g*sin(θ1 - 2θ2)
              - 2m2*sin(Δθ)*(θ̇2^2*l2 + θ̇1^2*l1*cos(Δθ)) ) / d1
              
    θ̈2 = ( 2*sin(Δθ)*(θ̇1^2*l1*(m1 + m2) + g*(m1 + m2)*cos(θ1)
              + θ̇2^2*l2*m2*cos(Δθ)) ) / d2
              
    return θ̈1, θ̈2
end

function residuals(x, θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt)
    θ1_next, θ2_next, θ̇1_next, θ̇2_next = x
    
    θ̈1_next, θ̈2_next = compute_accelerations(θ1_next, θ2_next, θ̇1_next, θ̇2_next)

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

function newton_raphson(θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt; max_iter=50, tol=1e-8)
    θ̈1_guess, θ̈2_guess = compute_accelerations(θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr)
    θ̇1_guess = θ̇1_curr + Δt * θ̈1_guess
    θ̇2_guess = θ̇2_curr + Δt * θ̈2_guess
    θ1_guess = θ1_curr + Δt * θ̇1_guess
    θ2_guess = θ2_curr + Δt * θ̇2_guess
    
    x = [θ1_guess, θ2_guess, θ̇1_guess, θ̇2_guess]
    
    for i in 1:max_iter
        r = residuals(x, θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt)
        
        if norm(r) < tol
            return x
        end
        
        J = jacobian(x, θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt)
        x -= J \ r
    end
    
    @warn "牛顿-拉夫逊方法未收敛，最大迭代次数已达: $max_iter"
    return x
end


for n in 1:nₜ-1
    θ1_curr = θ1[n]
    θ2_curr = θ2[n]
    θ̇1_curr = θ̇1[n]
    θ̇2_curr = θ̇2[n]

    x_next = newton_raphson(θ1_curr, θ2_curr, θ̇1_curr, θ̇2_curr, Δt)

    θ1[n+1], θ2[n+1], θ̇1[n+1], θ̇2[n+1] = x_next
end

fig = Figure(resolution=(1000, 600))
ax = Axis(fig[1, 1], xlabel="时间 (s)", ylabel="角度 (rad)", title="双摆角度随时间变化 (隐式欧拉法)")
lines!(ax, t, θ1, color=:blue, label="第一个摆角度 θ₁")
lines!(ax, t, θ2, color=:red, label="第二个摆角度 θ₂")
axislegend(ax)

fig
