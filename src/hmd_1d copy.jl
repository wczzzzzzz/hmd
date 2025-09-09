using  ApproxOperator, CairoMakie
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫∫q̇mṗqkpdxdt

using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

kᶜ = 100.0
m = 1.0
q̇₀ = 1.0
q₀ = 1.0
α = 1e12
β = 1e12
t = 0.0:0.1:8.0
𝜔 = (kᶜ/m)^0.5
𝑢(t) = q₀*cos(𝜔*t) + q̇₀/𝜔*sin(𝜔*t)
# 𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)

const to = TimerOutput()
gmsh.initialize()
# @timeit to "open msh file" gmsh.open("./msh/bar/bar_160.msh")
@timeit to "open msh file" gmsh.open("./msh/bar/bar_un_80.msh")
# @timeit to "open msh file" gmsh.open("./msh/bar/111.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nₚ = length(nodes)
k = zeros(nₚ,nₚ)
f = zeros(nₚ)

@timeit to "calculate ∫∫∇q∇pdxdt" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :m=>(x,y,z)->m, :kᶜ=>(x,y,z)->kᶜ)
    prescribe!(elements, :α=>(x,y,z)->α, :β=>(x,y,z)->β)
    # prescribe!(elements, :c=>(x,y,z)->c)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set𝝭!(elements)
    𝑎 = ∫∫q̇mṗqkpdxdt=>elements
    @timeit to "assemble" 𝑎(k)
end

𝑃₀ = m*q̇₀
f[1] -= 𝑃₀

α = 1e12
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵅ[1,1] += α
fᵅ[1] += α*q₀
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
kᵝ[2,2] += α
# kᵝ[nₚ,nₚ] += α

d = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
d = d[1:nₚ]
# e = d - 𝑢.(t)

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, t, d[[1,3:end...,2]], color = :blue)
lines!(ax, t, 𝑢, color = :black)
# lines!(nodes.x[[1,3:end...,2]], d[[1,3:end...,2]], color = :blue)
# lines!(nodes.x, d, color = :blue)
# lines!(t, e[[1,3:end...,2]], color = :red)


# e = d - 𝑥

# lines!(𝑡, e[[1,3:end...,2]], color = :red)


fig

# save("./fig/非均布弹簧小车/prog=1.05位移图/n=320.png",fig)
