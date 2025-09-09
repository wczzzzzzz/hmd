using  ApproxOperator

import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
using GLMakie
using WriteVTK
using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

# ps = MKLPardisoSolver()
# set_matrixtype!(ps,2)

ρA = 1.0
EA = 1.0
α = 1e15
L = 1.0
v₀ = 1.0

const to = TimerOutput()
gmsh.initialize()
@timeit to "open msh file" gmsh.open("./msh/b=2/Tri3反向20.msh")
# @timeit to "open msh file" gmsh.open("./msh/b=2/Tri6反向20.msh")
# @timeit to "open msh file" gmsh.open("./msh/b=2/Tri3非均布20.msh")
# @timeit to "open msh file" gmsh.open("./msh/b=2/Tri3非均布20.msh")

@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nₚ = length(nodes)
nₑ = length(elements)
k = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)

@timeit to "calculate ∫∫∇q∇pdxdt" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :EA=>(x,y,z)->EA, :ρA=>(x,y,z)->ρA)
    prescribe!(elements, :α=>(x,y,z)->α, :β=>(x,y,z)->β)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    𝑎 = ∫∫∇q∇pdxdt=>elements
    @timeit to "assemble" 𝑎(k)
end

@timeit to "calculate ∫vtdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"])
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"])
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"])
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"])
    prescribe!(elements_1,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_3,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_4,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_1,:t=>(x,y,z)->1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑓 = ∫vtdΓ=>elements_1
    𝑎ᵅ = ∫vgdΓ=>elements_1
    𝑎ᵅ = ∫vgdΓ=>elements_4
    𝑎ᵝ = ∫vgdΓ=>elements_3
    𝑎ᵝ = ∫vgdΓ=>elements_4
    @timeit to "assemble" 𝑓(f)
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
    @timeit to "assemble" 𝑎ᵝ(kᵝ,fᵝ)
end

@timeit to "solve" dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
d = dt[1:nₚ]
push!(nodes,:d=>d)

xs = [node.x for node in nodes]'
ys = [node.y for node in nodes]'
zs = [node.z for node in nodes]'
points = [xs; ys; zs]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements]

σ = zeros(nₑ)
for (j,p) in enumerate(elements)
    σ_ = 0.0
    𝑤_ = 0.0
    for ξ in p.𝓖
        B₁ = ξ[:∂𝝭∂x]
        ε = 0.0
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(p.𝓒)
            ε += B₁[i]*xᵢ.d
        end
        σ_ += EA*ε*𝑤
        𝑤_ += 𝑤
    end
    σ[j] = σ_/𝑤_
end

# vtk_grid("./vtk/Tri6_非均布_20"*string(nₚ), points, cells) do vtk
#     vtk["位移"] = d
#     vtk["应力"] = σ
# end
