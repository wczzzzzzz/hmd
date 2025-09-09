using ApproxOperator, XLSX, LinearAlgebra, LinearSolve, GLMakie
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁
import ApproxOperator.Hamilton: ∫qmpdΩ, ∫qkpdΩ

using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

ρA = 1.0
EA = 1.0
α = 1e9
β = 1e6
c = (EA/ρA)^0.5
𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)
function 𝑢(x,t)
    if x < t - 1
        return 2/π
    elseif x > t
        return 0.0
    else
        return (1-cos(π*(c*t - x)))/π
    end
end

const to = TimerOutput()
gmsh.initialize()
@timeit to "open msh file" gmsh.open("./msh/bar/bar_20.msh")

@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nₚ = length(nodes)
nₑ = length(elements)
k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
fᵗ = zeros(nₚ)
fᵍ = zeros(nₚ)
T = 4
Δt = 0.1
nₜ = Int(T/Δt)
d = zeros(nₚ,nₜ+1)
d̈ₙ₊₁ = zeros(nₚ)
ḋₙ = zeros(nₚ)
ḋₙ₊₁ = zeros(nₚ)

@timeit to "calculate qkpdΩ" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :EA=>(x,y,z)->EA, :ρA=>(x,y,z)->ρA)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set𝝭!(elements)
    𝑎 = ∫qkpdΩ=>elements
    b = ∫qmpdΩ=>elements
    @timeit to "assemble" 𝑎(k)
    @timeit to "assemble" b(m)
end

@timeit to "calculate ∫vtdΓ" begin
    @timeit to "get elements" elements_t = getElements(nodes, entities["Γᵗ"])
    @timeit to "get elements" elements_g = getElements(nodes, entities["Γᵍ"])
    prescribe!(elements_g,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    @timeit to "calculate shape functions" set𝝭!(elements_t)
    @timeit to "calculate shape functions" set𝝭!(elements_g)
    𝑎ᵍ = ∫vgdΓ=>elements_g
    @timeit to "assemble" 𝑎ᵍ(m,fᵍ)
end

for n in 1:nₜ
    fill!(fᵗ,0.0)
    t = (n+1)*Δt
    prescribe!(elements_t,:t=>(x,y,z)->-𝑇(t))
    𝑎ᵗ = ∫vtdΓ=>elements_t
    @timeit to "assemble" 𝑎ᵗ(fᵗ)
    d̈ₙ₊₁ .= m\(fᵗ+fᵍ - k*d[:,n])
    ḋₙ₊₁ .+= ḋₙ + Δt*d̈ₙ₊₁
    d[:,n+1] .= d[:,n] + Δt*ḋₙ₊₁
end
# push!(nodes,:d=>d[:,17])
# push!(nodes,:d=>d)

# elements_Ωᵍ = getElements(nodes, entities["Ω"], 10)
# prescribe!(elements_Ωᵍ,:u=>(x,y,z)->𝑢(x,y))
# prescribe!(elements_Ωᵍ,:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements_Ωᵍ,:∂u∂y=>(x,y,z)->∂u∂t(x,y))
# prescribe!(elements_Ωᵍ,:∂u∂z=>(x,y,z)->0.0)
# set∇𝝭!(elements_Ωᵍ)
# 𝐿₂ = log10.(L₂(elements_Ωᵍ))

ys = 0.0:4.0/(41-1):4.0
# ys = 0.0:Δt:T
XLSX.openxlsx("./excel/hmd_2d_Euler.xlsx", mode="w") do xf
Sheet = xf[1]
Sheet["A1"] = "x"
Sheet["B1"] = "t"
Sheet["C1"] = "d"
Sheet["D1"] = "误差Δ"
Sheet["E1"] = "log10(4/32)"
Sheet["F1"] = "𝐿₂"
row = 2
    for i in 1:nₚ
        x = nodes[i].x
        for j in 1:(nₜ+1)
            t = ys[j]
            z = d[i,j]
            Δ = z - 𝑢(x, t)
            Sheet["A$(row)"] = x
            Sheet["B$(row)"] = t
            Sheet["C$(row)"] = z
            Sheet["D$(row)"] = Δ
            Sheet["E$(row)"] = log10(4/32)
            Sheet["F$(row)"] = 𝐿₂
            row += 1
        end
    end
end


# fig = Figure()
# ax1 = Axis3(fig[1,1])

# xs = zeros(nₚ)
# ys = zeros(nₚ)
# ds = zeros(nₚ)
# us = zeros(nₚ)
# es = zeros(nₚ)

# for (i, node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     us[i] = 𝑢(x,y)
# end
# for (i,node) in enumerate(nodes)
#     xs[i] = node.x
#     ys[i] = node.y
#     ds[i] = node.d
#     es[i] = ds[i] - us[i]
# end
# face = zeros(nₑ,3)
# for (i,elm) in enumerate(elements)
#     face[i,:] .= [x.𝐼 for x in elm.𝓒]
# end

# meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)

# fig