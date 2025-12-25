using CairoMakie

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

dt = 0.01
t = collect(0.0:dt:10.0)
nₚ = length(t)
nₑ = nₚ-1

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

for i in 1:nₑ
    t₁ = t[i]
    t₂ = t[i+1]
    𝐿 = t₂ - t₁
    k[i,i] += l/𝐿 - g/3*𝐿
    k[i,i+1] += -l/𝐿 - g/6*𝐿
    k[i+1,i] += -l/𝐿 - g/6*𝐿
    k[i+1,i+1] += l/𝐿 - g/3*𝐿
end

𝑃₀ = l*θ̇₀
f[1] -= 𝑃₀

α = 1e12
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵅ[1,1] += α
fᵅ[1] += α*θ₀
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
# kᵝ[1,1] += α
kᵝ[nₚ,nₚ] += α

d = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
δd = d[nₚ+1:2*nₚ]
d = d[1:nₚ]

e = d - 𝑥
# lines!(t, e, color = :red)
lines!(t, d, color = :blue)
# lines!(t, δd, color = :red)

fig

# XLSX.openxlsx("./excel/摇摆问题.xlsx", mode="rw") do xf
#     Sheet = xf[4]
#     for i in 1:length(t)
#        Sheet["A$(i)"] = t[i]
#        Sheet["B$(i)"] = d[i]
#        Sheet["C$(i)"] = 𝑥[i]
#        Sheet["D$(i)"] = e[i]
#     end
# end