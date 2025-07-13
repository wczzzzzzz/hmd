using GLMakie

k = 100
m = 1.0
q̇₀ = 1.0
q₀ = 1.0

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "T", ylabel = "u")
𝜔 = (k/m)^0.5
𝑢(t) = q₀*cos(𝜔*t) + q̇₀/𝜔*sin(𝜔*t)

T = 8.0
Δt = 0.005
t = 0.0:Δt:T
nₜ = Int(T/Δt)

q = zeros(nₜ+1)
q[1] = q₀
q̇ = zeros(nₜ+1)
q̇[1] = q̇₀

for i in 1:nₜ
    q̈[i] = -k/m * q[i]
    q̇[i+1] = q̇[i] + Δt * q̈[i]
    q[i+1] = q[i] + Δt * q̇[i]
end

invisible_line = lines!(ax, [0, 0], [0, 0], color = :white, label="Δt=0.005", visible=false)
blue_line = lines!(ax, t[1:2], q[1:2], color = :blue, label="Explicit Euler")
black_line = lines!(ax, t[1:2], 𝑢.(t[1:2]), color = :black, label="Exact Solution")
leg = Legend(fig[1, 2], [blue_line, black_line, invisible_line], ["Explicit Euler", "Exact Solution", "Δt=0.005"], position=(0.95, 0.95))

nᵢ = 10
record(fig, "./fig/一维/Explicit_1d.gif", 1:nₜ,framerate = 50) do i
    t = 0.0:Δt:i*Δt
    lines!(ax, t, q[1:i+1], color = :blue)
    t = 0.0:Δt/nᵢ:i*Δt
    lines!(ax, t, 𝑢.(t), color = :black)
end

# e = q - 𝑥
# lines!(t, e, color = :red)

# save("./fig/一维/Explicit_1d_e.png",fig)

xlims!(ax, 0, 8)
fig