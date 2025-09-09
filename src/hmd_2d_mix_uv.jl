using ApproxOperator, XLSX, LinearAlgebra, LinearSolve, GLMakie
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, stabilization_bar_LSG

using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

ρA = 1.0
EA = 1.0
α = 1e6
𝑇(t) = t > 1.0 ? 0.0 : - sin(π*t)
function 𝑢(x,t)
    if x < t - 1
        return 2/π
    elseif x > t
        return 0.0
    else
        return (1-cos(π*(t - x)))/π
    end
end
function ∂u∂t(x, t)
    if x < t - 1 || x > t
        return 0.0
    else
        return sin(π * (t - x))
    end
end
function ∂u∂x(x, t)
    if x < t - 1
        return 0.0
    elseif x > t
        return 0.0
    else
        return -sin(π*(t - x))
    end
end

const to = TimerOutput()
gmsh.initialize()
# @timeit to "open msh file" gmsh.open("./msh/square/Tri6_16.msh")
@timeit to "open msh file" gmsh.open("./msh/Non-uniform/Tri6/8.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()
# @timeit to "open msh file" gmsh.open("./msh/square/Tri3_16.msh")
@timeit to "open msh file" gmsh.open("./msh/Non-uniform/Tri3/8.msh")
@timeit to "get entities" entities_p = getPhysicalGroups()
@timeit to "get nodes" nodes_p = get𝑿ᵢ()

nᵤ = length(nodes)
nₚ = length(nodes_p)
# nₑ = length(elements)
kᵤᵤ = zeros(nᵤ,nᵤ)
kᵅ = zeros(nᵤ,nᵤ)
fᵅ = zeros(nᵤ)
kᵤₚ = zeros(nᵤ,nₚ)
fₚ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

@timeit to "calculate ∫∫∇q∇pdxdt" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"], 8)
    @timeit to "get elements" elements_p = getElements(nodes_p, entities_p["Ω"], 8)
    prescribe!(elements, :EA=>(x,y,z)->EA, :ρA=>(x,y,z)->ρA)
    prescribe!(elements_p, :EA=>(x,y,z)->EA, :ρA=>(x,y,z)->ρA)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set∇𝝭!(elements_p)
    𝑎ᵘ = ∫∫∇q∇pdxdt=>elements
    𝑎ᵘᵖ = ∫∫∇q∇pdxdt=>(elements,elements_p)
    @timeit to "assemble" 𝑎ᵘ(kᵤᵤ)
    @timeit to "assemble" 𝑎ᵘᵖ(kᵤₚ)
end

@timeit to "calculate ∫vtdΓ, ∫vgdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"], 8)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"], 8)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"], 8)
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"], 8)
    @timeit to "get elements" elements_1p = getElements(nodes_p, entities_p["Γ¹"], 8)
    @timeit to "get elements" elements_2p = getElements(nodes_p, entities_p["Γ²"], 8)
    @timeit to "get elements" elements_3p = getElements(nodes_p, entities_p["Γ³"], 8)
    @timeit to "get elements" elements_4p = getElements(nodes_p, entities_p["Γ⁴"], 8)
    prescribe!(elements_1,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_2,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_3,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_2p,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_3p,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_4p,:t=>(x,y,z)->-𝑇(y))
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    @timeit to "calculate shape functions" set𝝭!(elements_1p)
    @timeit to "calculate shape functions" set𝝭!(elements_2p)
    @timeit to "calculate shape functions" set𝝭!(elements_3p)
    @timeit to "calculate shape functions" set𝝭!(elements_4p)
    𝑓ᵖ = ∫vtdΓ=>elements_4p
    𝑎ᵅ = ∫vgdΓ=>elements_1
    𝑎ᵅ = ∫vgdΓ=>elements_2
    𝑎ᵝ = ∫vgdΓ=>elements_2p
    𝑎ᵝ = ∫vgdΓ=>elements_3p
    @timeit to "assemble" 𝑓ᵖ(fₚ)
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
    @timeit to "assemble" 𝑎ᵝ(kᵝ,fᵝ)
end

@timeit to "solve" dt = [kᵤᵤ+kᵅ -kᵤₚ;-kᵤₚ' kᵝ]\[fᵅ;-fₚ+fᵝ]
d = dt[1:nₚ]
# δd = dt[nₚ+1:end]
push!(nodes, :d=>d)
# push!(nodes,:δd=>δd)

# elements_Ωᵍ = getElements(nodes, entities["Ω"], 10)
# prescribe!(elements_Ωᵍ,:u=>(x,y,z)->𝑢(x,y))
# prescribe!(elements_Ωᵍ,:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements_Ωᵍ,:∂u∂y=>(x,y,z)->∂u∂t(x,y))
# prescribe!(elements_Ωᵍ,:∂u∂z=>(x,y,z)->0.0)
# set∇𝝭!(elements_Ωᵍ)
# 𝐿₂ = log10.(L₂(elements_Ωᵍ))
# println(𝐿₂)

# gmsh.finalize()

fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nᵤ)
ys = zeros(nᵤ)
ds = zeros(nᵤ)
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    ds[i] = node.d
end
xp = zeros(nₚ)
yp = zeros(nₚ)
δds = zeros(nₚ)
for (i,node) in enumerate(nodes_p)
    xp[i] = node.x
    yp[i] = node.y
    # δds[i] = node.δd
end
face = zeros(nₑ,6)
for (i,elm) in enumerate(elements)
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax,xs,ys,zs,color=zs,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax2,xp,yp,δds,color=δds,markersize = 0.1)
fig

# save("./fig/616测试/非均布_32.png",fig)

# index = [4,8,16,32]
# XLSX.openxlsx("./excel/hmd_2d_mix_uv.xlsx", mode="rw") do xf
#     Sheet = xf[1]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = log10(4/ndiv)
#     Sheet["B"*string(ind)] = 𝐻₁
#     Sheet["C"*string(ind)] = 𝐿₂
# end
 