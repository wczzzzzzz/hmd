using GLMakie
using LinearAlgebra


m = 1.0
k = 100.0
b = 0.5
Γ = b/m 
ω₀ = sqrt(k/m)
Δt = 0.01
T = 8.0
t = 0.0:Δt:T
nₜ = Int(T/Δt)

x₀ = 1.0 
v₀ = 1.0

x = zeros(nₜ+1)
v = zeros(nₜ+1)
a = zeros(nₜ+1)
x[1] = x₀
v[1] = v₀
    
for n in 1:nₜ
    a[n] = -Γ*v[n]-ω₀^2*x[n]
    v[n+1] = v[n] + Δt * a[n]
    x[n+1] = x[n] + Δt * v[n+1]
end

# for n in 1:nₜ
#     d = 1 + Δt * Γ + (Δt * ω₀)^2
#     v[n+1] = (v[n] - Δt * ω₀^2 * x[n]) / d
#     x[n+1] = x[n] + Δt * v[n+1]
# end

fig = Figure()
ax = Axis(fig[1, 1])
𝑡 = 0.0:0.01:8.0
ω₁ = sqrt(ω₀^2-(Γ*Γ/4))
# 𝑥 = e^(-Γ*t/2)(x₀*cos.(ω₁*t)+(v₀+Γ*x₀/2)*sin.(ω₁*t)/ω₁)

lines!(ax, t, x, color = :blue)
# lines!(ax, t, 𝑥, color = :black)


fig

