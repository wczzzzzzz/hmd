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
@timeit to "open msh file" gmsh.open("./msh/bar/bar_21.msh")
# @timeit to "open msh file" gmsh.open("./msh/b=2/Tri6反向20.msh")
# @timeit to "open msh file" gmsh.open("./msh/b=2/Tri3非均布20.msh")
# @timeit to "open msh file" gmsh.open("./msh/b=2/Tri3非均布20.msh")

@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nₚ = length(nodes)
k = zeros(nₚ,nₚ)
m = zeros(nₚ,nₚ)
fᵗ = zeros(nₚ)
fᵍ = zeros(nₚ)
T = 2
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
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"])
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"])
    prescribe!(elements_1,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_1,:t=>(x,y,z)->1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    𝑎ᵅ = ∫vgdΓ=>elements_1
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
end

for n in 1:nₜ
    fill!(fᵗ,0.0)
    prescribe!(elements_t,:t=>(x,y,z)->-1.0)
    𝑎ᵗ = ∫vtdΓ=>elements_t
    @timeit to "assemble" 𝑎ᵗ(fᵗ)
    d̈ₙ₊₁ .= m\(fᵗ+fᵍ - k*d[:,n])
    ḋₙ₊₁ .+= ḋₙ + Δt*d̈ₙ₊₁
    d[:,n+1] .= d[:,n] + Δt*ḋₙ₊₁
end

push!(nodes,:d=>d[:,11])

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
# vtk_grid("./vtk/Tri6_非均布_"*string(ndiv)*"_"*string(nₚ), points, cells) do vtk
#     vtk["位移"] = d
#     vtk["应力"] = σ
# end
