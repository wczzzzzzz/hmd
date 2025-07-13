using   GLMakie

k = 100
m = 1.0
q̇₀ = 1.0
q₀ = 1.0

fig = Figure()
# ax = Axis(fig[1, 1])
ax = Axis(fig[1, 1], xlabel = "T", ylabel = "x")
𝑡 = 0.0:0.005:8.0
𝜔 = (k/m)^0.5
𝑢(t) = q₀*cos(𝜔*t) + q̇₀/𝜔*sin(𝜔*t)
# 𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(ax, 𝑡, 𝑢, color = :black)

T = 8.0
Δt = 0.005
t = 0.0:Δt:T
nₜ = Int(T/Δt)

q = zeros(nₜ+1)
q̇ = zeros(nₜ+1)
q̈ = zeros(nₜ+1)

q[1] = q₀
q̇[1] = q̇₀

for n in 1:nₜ
    q̈[n] = -k/m * q[n]  
    q̇[n+1] = q̇[n] + Δt * q̈[n]  
    q[n+1] = q[n] + Δt * q̇[n]  
end

invisible_line = lines!(ax, [0, 0], [0, 0], color = :white, label="Δt=0.005", visible=false)
blue_line = lines!(ax, t[1:2], q[1:2], color = :blue, label="Explicit Euler")
black_line = lines!(ax, t[1:2], 𝑢.(t[1:2]), color = :black, label="Exact Solution")
# red_line = lines!(ax, t[1:2], e[1:2], color = :red, label="error")
leg = Legend(fig[1, 2], [blue_line, black_line, invisible_line], ["Explicit Euler", "Exact Solution", "Δt=0.005"], position=(0.95, 0.95))
# leg = Legend(fig[1, 2], [red_line, invisible_line], ["error", "Δt=0.005"], position=(0.95, 0.95))

lines!(ax, t, q, color = :blue)

# e = q - 𝑢
# lines!(t, e, color = :red)
xlims!(ax, 0, 8)

fig

# save("./fig/一维/Explicit_1d.png",fig)
# save("./fig/一维/Explicit_1d_e.png",fig)


