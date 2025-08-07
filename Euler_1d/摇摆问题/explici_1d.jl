using GLMakie

g = 9.81  
l = 1.0   
θ₀ = 0.1
θ̇₀ = 0.0 

fig = Figure()
ax = Axis(fig[1, 1])
𝑡 = 0.0:0.01:10.0
𝜔 = (g/l)^0.5
# 𝑢(t) = θ₀*cos(𝜔*𝑡) + (θ̇₀/𝜔)*sin(𝜔*𝑡)
# 𝑥 = 𝑢.(𝑡)
𝑥 = θ₀.*cos.(𝜔.*𝑡) + θ̇₀/𝜔.*sin.(𝜔.*𝑡)
lines!(ax, 𝑡, 𝑥, color = :black)

T = 10.0 
Δt = 0.01 
t = 0.0:Δt:T
nₜ = Int(T/Δt)

θ = zeros(nₜ+1)
θ̇  = zeros(nₜ+1)
θ̈  = zeros(nₜ+1)

θ[1] = θ₀
θ̇[1] = θ̇₀

for n in 1:nₜ
    θ̈[n] = -g/l * θ[n]  
    θ̇[n+1] = θ̇[n] + Δt * θ̈[n]
    θ[n+1] = θ[n] + Δt * θ̇[n]
end

lines!(ax, t, θ, color = :blue)
e = θ - 𝑥
# xlims!(ax, 0, 10)
lines!(t, e, color = :red)

fig

# XLSX.openxlsx("./excel/摇摆问题.xlsx", mode="rw") do xf
#     Sheet = xf[1]
#     for i in 1:length(t)
#        Sheet["A$(i)"] = t[i]
#        Sheet["B$(i)"] = θ[i]
#        Sheet["C$(i)"] = 𝑥[i]
#        Sheet["D$(i)"] = e[i]
#     end
# end