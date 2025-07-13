using ApproxOperator, XLSX, LinearAlgebra, LinearSolve, GLMakie
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, stabilization_bar_LSG, stabilization_bar_LSG_Γ

using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

ρA = 1.0
EA = 1.0
α = 1e6
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
function ∂u∂t(x, t)
    if x < t - 1 || x > t
        return 0.0
    else
        return sin(π * (c*t - x))
    end
end
function ∂u∂x(x, t)
    if x < t - 1
        return 0.0
    elseif x > t
        return 0.0
    else
        return -sin(π*(c*t - x))
    end
end

const to = TimerOutput()
gmsh.initialize()
# @timeit to "open msh file" gmsh.open("./msh/Non-uniform/Tri6/16.msh")
# @timeit to "open msh file" gmsh.open("./msh/Non-uniform/Tri3/4.msh")
# @timeit to "open msh file" gmsh.open("./msh/Non-uniform/拉伸压缩C=1.0/2.0_4.msh")
@timeit to "open msh file" gmsh.open("./msh/square/square_16.msh")
# @timeit to "open msh file" gmsh.open("./msh/square/Tri6_4")

@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nₚ = length(nodes)
nₑ = length(elements)
k = zeros(nₚ,nₚ)
kˢ = zeros(nₚ,nₚ)
f = zeros(nₚ)
kᵅ = zeros(nₚ,nₚ)
fᵅ = zeros(nₚ)
kᵝ = zeros(nₚ,nₚ)
fᵝ = zeros(nₚ)
kᵞ = zeros(nₚ,nₚ)
kᵗ = zeros(nₚ,nₚ)

@timeit to "calculate ∫∫∇q∇pdxdt" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :EA=>(x,y,z)->EA, :ρA=>(x,y,z)->ρA)
    prescribe!(elements, :α=>(x,y,z)->α, :β=>(x,y,z)->β)
    prescribe!(elements, :c=>(x,y,z)->c)
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
    prescribe!(elements_2,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_3,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_4,:t=>(x,y,z)->-𝑇(y))
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑓 = ∫vtdΓ=>elements_4
    𝑎ᵅ = ∫vgdΓ=>elements_1
    𝑎ᵅ = ∫vgdΓ=>elements_2
    𝑎ᵝ = ∫vgdΓ=>elements_2
    𝑎ᵝ = ∫vgdΓ=>elements_3
    @timeit to "assemble" 𝑓(f)
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
    @timeit to "assemble" 𝑎ᵝ(kᵝ,fᵝ)
end

@timeit to "solve" dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
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

fig = Figure()
ax1 = Axis3(fig[1,1])
# ax2 = Axis3(fig[1,2])

xs = zeros(nₚ)
ys = zeros(nₚ)
ds = zeros(nₚ)
# δds = zeros(nₚ)
us = zeros(nₚ)
# qs = zeros(nₚ)
# as = zeros(nₚ)
es = zeros(nₚ)

for (i, node) in enumerate(nodes)
    x = node.x
    y = node.y
    us[i] = 𝑢(x,y)
    # qs[i] = ∂u∂t(x,y)
    # as[i] = ∂²u∂t²(x,y)
end
for (i,node) in enumerate(nodes)
    xs[i] = node.x
    ys[i] = node.y
    ds[i] = node.d
    # δds[i] = node.δd
    es[i] = ds[i] - us[i]
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements)
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# # mesh!(ax,xs,ys,zs,face,color=ds)
# # meshscatter!(ax1,xs,ys,us,color=us,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# # meshscatter!(ax1,xs,ys,es,color=es,markersize = 0.06)
# meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.06)
fig

# save("./fig/hmd_2d/test_x=20/t=98.png",fig)
# save("./fig/72测试/Tri6_非均布_LSG_32.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri3/三维图/t=25.png",fig)
# save("./fig/hmd_2d/锁三边x=20/Tri6/均布/t=25.png",fig)
# save("./fig/hmd_2d/局部加密C=0.2/T6_c=0.05.png",fig)
# save("./fig/hmd_2d/Tri3/非均布/n=80.png",fig)

# for i in 1:nₚ
#     x = nodes.x[i]
#     y = nodes.y[i]
#     d₁ = d[i]
#     # Δ = d[i] - 𝑢(x,y)
#         index = [10,20,40,80]
#         XLSX.openxlsx("./excel/square.xlsx", mode="rw") do xf
#         Sheet = xf[3]
#         ind = findfirst(n->n==ndiv,index)+1
#         # Sheet["A"*string(ind)] = x
#         # Sheet["B"*string(ind)] = y
#         # Sheet["C"*string(ind)] = d₁
#         # Sheet["D"*string(ind)] = Δ
#         Sheet["E"*string(ind)] = 𝐿₂
#         Sheet["F"*string(ind)] = log10(4/ndiv)
#     end
# end

# index = [4,8,16,32]
# # index = [0.4,0.3,0.2,0.1]
# # index = [0,1,2,3]
# XLSX.openxlsx("./excel/hmd_2d.xlsx", mode="rw") do xf
#     Sheet = xf[3]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = log10(4/ndiv)
#     # Sheet["A"*string(ind)] = log10(nₚ)
#     # Sheet["B"*string(ind)] = 𝐻₁
#     Sheet["C"*string(ind)] = 𝐿₂
# end

# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = node.d
#     # points[3,i] = us[i]*4
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/hmd_2d/error/non_uniform_Tri3_"*string(ndiv)*".vtu",points,cells) do vtk
# vtk_grid("./vtk/hmd_2d/exact_d_"*string(ndiv)*".vtu",points,cells) do vtk
#     # vtk["d"] = [node.d for node in nodes]
#     vtk["精确解"] = us
# end

# fₓ,fₜ,fₓₓ,fₜₜ = truncation_error(elements["Ω"],nₚ)
# println(fₓ)
# println(fₜ)
# println(fₛ)

# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# vtk_grid("./vtk/hmd_2d/error/uniform_Tri3_"*string(ndiv), points, cells) do vtk
#     # vtk["fₓ"] = fₓ
#     # vtk["fₜ"] = fₜ
#     # vtk["fₓₓ"] = fₓₓ
#     # vtk["fₜₜ"] = fₜₜ
#     # vtk["fₓₓ/fₜₜ"] = fₓₓ./fₜₜ
#     vtk["误差"] = es
# end

