using ApproxOperator, XLSX, LinearAlgebra, LinearSolve, GLMakie
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, stabilization_bar_LSG

using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

α = 1e6
β = 1e6
# ρA = 1.0*400.0/100.0
ρA = 1.0
EA = 1.0
a = 1.0
l = 4.0
c = (EA/ρA)^0.5
φ(x) = sin(π*x/l)
𝑢(x,t) = cos.(π.*a.*c*t/l).*sin.(π.*x/l)
∂u∂t(x,t) = (-π.*a.*c/l)*sin.(π.*a.*c*t/l).*sin.(π.*x/l)
∂u∂x(x,t) = (π./l)*cos.(π.*a.*c*t/l).*cos.(π.*x/l)

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
k = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
kˢ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
f = zeros(nₚ + nₗ + nₑ)
kᵅ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
fᵅ = zeros(nₚ + nₗ + nₑ)
kᵝ = spzeros(nₚ + nₗ + nₑ, nₚ + nₗ + nₑ)
fᵝ = zeros(nₚ + nₗ + nₑ)

@timeit to "calculate ∫∫∇q∇pdxdt" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    @timeit to "get elements" elements_t, nodes_t, edges = Tri3toTriHermite(elements, nodes)
    prescribe!(elements_t, :EA=>(x,y,z)->EA, :ρA=>(x,y,z)->ρA)
    @timeit to "calculate shape functions" set𝝭!(elements_t)
    @timeit to "calculate shape functions" set∇𝝭!(elements_t)
    𝑎 = ∫∫∇q∇pdxdt=>elements_t
    @timeit to "assemble" 𝑎(k)
end

@timeit to "calculate ∫vtdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"])
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"])
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"])
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"])
    @timeit to "get elements" elements_1t = Seg2toSegHermite(elements_1, nodes_t, edges)
    @timeit to "get elements" elements_2t = Seg2toSegHermite(elements_2, nodes_t, edges)
    @timeit to "get elements" elements_3t = Seg2toSegHermite(elements_3, nodes_t, edges)
    @timeit to "get elements" elements_4t = Seg2toSegHermite(elements_4, nodes_t, edges)
    prescribe!(elements_1t,:t=>(x,y,z)->0.0)
    prescribe!(elements_1t,:g=>(x,y,z)->φ(x), :α=>(x,y,z)->α)
    prescribe!(elements_2t,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_3t,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_4t,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    @timeit to "calculate shape functions" set𝝭!(elements_1t)
    @timeit to "calculate shape functions" set𝝭!(elements_2t)
    @timeit to "calculate shape functions" set𝝭!(elements_3t)
    @timeit to "calculate shape functions" set𝝭!(elements_4t)
    𝑓 = ∫vtdΓ=>elements_1t
    𝑎ᵅ = ∫vgdΓ=>elements_1t
    𝑎ᵅ = ∫vgdΓ=>elements_2t
    𝑎ᵅ = ∫vgdΓ=>elements_4t
    𝑎ᵝ = ∫vgdΓ=>elements_2t
    𝑎ᵝ = ∫vgdΓ=>elements_3t
    𝑎ᵝ = ∫vgdΓ=>elements_4t
    @timeit to "assemble" 𝑓(f)
    @timeit to "assemble" 𝑎ᵅ(kᵅ,fᵅ)
    @timeit to "assemble" 𝑎ᵝ(kᵝ,fᵝ)
end

@timeit to "solve" dt = [k+kᵅ -k;-k kᵝ]\[fᵅ;-f+fᵝ]
# @timeit to "solve" dt = sparse([k+kᵅ -k;-k kᵝ])\[fᵅ;-f+fᵝ]
d = dt[1:nₚ]
# d = dt[1:nₚ+nₗ+nₑ]
push!(nodes, :d=>d)


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

xs = zeros(nₚ)
ys = zeros(nₚ)
zs = zeros(nₚ)
ds = zeros(nₚ)
# δds = zeros(nₚ)
# es = zeros(nₚ)
us = zeros(nₚ)
# qs = zeros(nₚ)
# as = zeros(nₚ)
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
    # es[i] = ds[i] - us[i]
end
face = zeros(nₑ,3)
for (i,elm) in enumerate(elements)
    face[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# mesh!(ax,xs,ys,face,color=zs)
# meshscatter!(ax1,xs,ys,es,color=es,markersize = 0.1)
meshscatter!(ax1,xs,ys,ds,color=ds,markersize = 0.06)
# meshscatter!(ax1,xs,ys,us,color=us,markersize = 0.1)
# meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.06)
fig

# # save("./fig/测试/74测试/Tri6_非均布_LSG_c^2_4.png",fig)

# # save("./fig/连续解/锁时间末端Tri_6非均布/t=19.png",fig)
# # save("./fig/连续解/锁时间末端Tri_6均布/t=25.png",fig)
# # save("./fig/连续解/mix_Tri_6均布/t=25.png",fig)
# # save("./fig/连续解/mix_Tri_6非均布/n=41.png",fig)

# # index = [4,8,16,32]
# # # index = [5,10,20,40]
# # XLSX.openxlsx("./excel/hmd_Continuous(2).xlsx", mode="rw") do xf
# #     Sheet = xf[5]
# #     ind = findfirst(n->n==ndiv,index)+1
# #     Sheet["A"*string(ind)] = log10(4/ndiv)
# #     # Sheet["A"*string(ind)] = log10(nₚ)
# #     # Sheet["B"*string(ind)] = 𝐻₁
# #     Sheet["C"*string(ind)] = 𝐿₂
# # end

# # points = zeros(3,nₚ)
# # for (i,node) in enumerate(nodes)
# #     points[1,i] = node.x
# #     points[2,i] = node.y
# #     points[3,i] = node.d
# # end
# # cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/hmd_Continuous/uniform_"*string(ndiv)*".vtu",points,cells) do vtk
# #     vtk["d"] = [node.d for node in nodes]
# # end

# # xs = [node.x for node in nodes]'
# # ys = [node.y for node in nodes]'
# # zs = [node.z for node in nodes]'
# # points = [xs; ys; zs]
# # cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/hmd_Continuous/error_uniform_"*string(ndiv), points, cells) do vtk
# #     vtk["误差"] = es
# #     # vtk["二阶导"] = as
# # end
