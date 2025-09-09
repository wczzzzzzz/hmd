using CairoMakie

k = 100
m = 1.0
q̇₀ = 1.0
q₀ = 1.0

fig = Figure(resolution = (700, 510))

ax = Axis(fig[1, 1], xlabel = "T", ylabel = "x", width = 600, height = 400, xticklabelsize = 18, yticklabelsize = 18)  

𝑡 = 0.0:0.005:8.0
𝜔 = (k / m)^0.5

𝑢(t) = q₀ * cos(𝜔 * t) + q̇₀ / 𝜔 * sin(𝜔 * t)

lines!(ax, 𝑡, 𝑢, color = :black)

T = 8.0
Δt = 0.005
t = 0.0:Δt:T
nₜ = Int(T / Δt)
q = zeros(nₜ + 1)
q̇ = zeros(nₜ + 1)
q̈ = zeros(nₜ + 1)

q[1] = q₀
q̇[1] = q̇₀

# for n in 1:nₜ
#     q̈[n] = -k / m * q[n]
#     q̇[n + 1] = q̇[n] + Δt * q̈[n]
#     q[n + 1] = q[n] + Δt * q̇[n + 1]
# end

# for n in 1:nₜ
#     q̈[n] = -k/m * q[n]  
#     q̇[n+1] = q̇[n] + Δt * q̈[n]  
#     q[n+1] = q[n] + Δt * q̇[n]  
# end

for n in 1:nₜ
    q[n+1] = (m*q[n])/(m + k*Δt^2) + (m*Δt*q̇[n])/(m + k*Δt^2)
    q̈[n+1] = -k/m * q[n+1]  
    q̇[n+1] = q̇[n] + Δt * q̈[n+1]  
end


invisible_line = lines!(ax, [0, 0], [0, 0], color = :white, label="Δt=0.005", visible=false)
blue_line = lines!(ax, t, q, color = :blue, label="Implicit Euler", linewidth = 3.0)
black_line = lines!(ax, t, 𝑢.(t), color = :black, label="Exact Solution")
ylims!(ax, -1.25, 1)
# ylims!(ax, -8, 7)

Legend(
    fig,  
    [blue_line, black_line, invisible_line], 
    ["Implicit Euler", "Exact Solution", "Δt=0.005"],
    halign = :left, 
    valign = :bottom,
    orientation = :horizontal,
    margin = (88, 10, 80, 10),
    tellheight = false,
    tellwidth = false,
    position = (0.05, -0.05)
)

xlims!(ax, 0, 8)

fig
# save("./fig/一维/semi_word_1d.png",fig)
# save("./fig/一维/explicit_word_1d.png",fig)
# save("./fig/一维/implicit_word_1d.png",fig)
