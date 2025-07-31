using  ApproxOperator, CairoMakie
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, ∫∫q̇mṗqkpdxdt

using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

kᶜ = 100
m = 1.0
q̇₀ = 1.0
q₀ = 1.0
α = 1e12
𝑡 = 0.0:0.005:8.0
𝜔 = (kᶜ/m)^0.5
𝑢(t) = q₀*cos(𝜔*t) + q̇₀/𝜔*sin(𝜔*t)
# 𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)

const to = TimerOutput()
gmsh.initialize()
@timeit to "open msh file" gmsh.open("./msh/bar/bar_160.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

@timeit to "calculate ∫∫∇q∇pdxdt" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :m=>(x,y,z)->m, :kᶜ=>(x,y,z)->kᶜ)
    prescribe!(elements, :α=>(x,y,z)->α, :β=>(x,y,z)->β)
    prescribe!(elements, :c=>(x,y,z)->c)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set𝝭!(elements)
    𝑎 = ∫∫q̇mṗqkpdxdt=>elements
    @timeit to "assemble" 𝑎(k)

    @timeit to "get elements" elements_t = getElements(nodes, entities["Γᵗ"])
    @timeit to "get elements" elements_g = getElements(nodes, entities["Γᵍ"])
    prescribe!(elements_t,:g=>(x,y,z)->1.0, :α=>(x,y,z)->α)
    prescribe!(elements_g,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_t,:t=>(x,y,z)->m*q̇₀)
    @timeit to "calculate shape functions" set𝝭!(elements_t)
    @timeit to "calculate shape functions" set𝝭!(elements_g)
    𝑓 = ∫vtdΓ=>elements_t
    𝑎ᵅ = ∫vgdΓ=>elements_t
    # 𝑎ᵅ = ∫vgdΓ=>elements_g
    # 𝑎ᵝ = ∫vgdΓ=>elements_t
    𝑎ᵝ = ∫vgdΓ=>elements_g
    @timeit to "assemble" 𝑓(f)
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
    @timeit to "assemble" 𝑎ᵝ(kᵝ,fᵝ)
end

@timeit to "solve" dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
d = dt[1:nₚ]
# δd = dt[nₚ+1:end]
push!(nodes, :d=>d)



# include("import_hmd.jl")

# ndiv= 160
# elements,nodes = import_hmd_bar("./msh/bar/bar_"*string(ndiv)*".msh")
# nₚ = length(nodes)
# nₑ = length(elements["Ω"])

# set𝝭!(elements["Ω"])
# set∇𝝭!(elements["Ω"])
# set𝝭!(elements["Γᵍ"])

# kᶜ = 100
# m = 1.0
# q̇₀ = 1.0
# q₀ = 1.0
# prescribe!(elements["Ω"],:m=>(x,y,z)->m)
# prescribe!(elements["Ω"],:kᶜ=>(x,y,z)->kᶜ)

# fig = Figure()
# Axis(fig[1, 1])
# 𝑡 = 0.0:0.005:8.0
# 𝜔 = (kᶜ/m)^0.5
# # 𝑢(t) = q₀*cos(𝜔*t) + q̇₀/𝜔*sin(𝜔*t)
# 𝑥 = q₀.*cos.(𝜔.*𝑡) + q̇₀/𝜔.*sin.(𝜔.*𝑡)
# # lines!(𝑡, 𝑥, color = :black)

# k = zeros(nₚ,nₚ)
# f = zeros(nₚ)

# 𝑎 = ∫∫q̇mṗqkpdxdt=>elements["Ω"]

# 𝑎(k)

# 𝑃₀ = m*q̇₀
# f[1] -= 𝑃₀

# α = 1e12
# kᵅ = zeros(nₚ,nₚ)
# fᵅ = zeros(nₚ)
# kᵅ[1,1] += α
# fᵅ[1] += α*q₀
# kᵝ = zeros(nₚ,nₚ)
# fᵝ = zeros(nₚ)
# # kᵝ[1,1] += α
# kᵝ[nₚ,nₚ] += α

# d = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# d = d[1:nₚ]


# lines!(nodes.x[[1,3:end...,2]], d[[1,3:end...,2]], color = :blue)
# # lines!(nodes.x, d, color = :blue)

# # e = d - 𝑢.(t)
# # e = d - 𝑥

# # lines!(𝑡, e[[1,3:end...,2]], color = :red)


# fig

# save("./fig/一维/hmd_1d.png",fig)
