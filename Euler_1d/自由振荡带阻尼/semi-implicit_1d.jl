using GLMakie

g = 9.81 
l1 = 1.0
l2 = 1.0
m1 = 1.0
m2 = 1.0

θ1₀ = 0.5
θ2₀ = 0.3
θ̇1₀ = 0.0
θ̇2₀ = 0.0

fig = Figure()
# ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

T = 20.0  
Δt = 0.01
t = 0.0:Δt:T
# nₜ = length(t)
nₜ = Int(T/Δt)

θ1 = zeros(nₜ+1)
θ2 = zeros(nₜ+1)
θ̇1 = zeros(nₜ+1)
θ̇2 = zeros(nₜ+1)
θ̈1 = zeros(nₜ+1)
θ̈2 = zeros(nₜ+1)

θ1[1] = θ1₀
θ2[1] = θ2₀
θ̇1[1] = θ̇1₀
θ̇2[1] = θ̇2₀

for n in 1:nₜ
    Δθ = θ1[n] - θ2[n]
    d1 = (m1 + m2)*l1 - m2*l1*cos(Δθ)^2
    d2 = (l2/l1)*d1
    θ̈1[n] = ( -g*(2m1 + m2)*sin(θ1[n]) - m2*g*sin(θ1[n] - 2θ2[n])
              - 2m2*sin(Δθ)*(θ̇2[n]^2*l2 + θ̇1[n]^2*l1*cos(Δθ)) ) / d1
    θ̈2[n] = ( 2*sin(Δθ)*(θ̇1[n]^2*l1*(m1 + m2) + g*(m1 + m2)*cos(θ1[n])
              + θ̇2[n]^2*l2*m2*cos(Δθ)) ) / d2
    θ̇1[n+1] = θ̇1[n] + Δt * θ̈1[n]
    θ̇2[n+1] = θ̇2[n] + Δt * θ̈2[n]
    θ1[n+1] = θ1[n] + Δt * θ̇1[n+1]
    θ2[n+1] = θ2[n] + Δt * θ̇2[n+1]
end

ax1 = Axis(fig[1, 1], xlabel = "时间 (s)", ylabel = "角度 (rad)", title = "双摆角度随时间变化")
lines!(ax1, t, θ1, color = :blue, label = "第一个摆角度 θ₁")
lines!(ax1, t, θ2, color = :red, label = "第二个摆角度 θ₂")
axislegend(ax1)

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

