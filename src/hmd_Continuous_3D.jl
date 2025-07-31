using ApproxOperator, XLSX, LinearAlgebra, LinearSolve, GLMakie
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫vtdΓ, ∫vgdΓ, ∫vbdΩ, L₂, ∫∫∇v∇udxdy, H₁
import ApproxOperator.Hamilton: ∫∫∇q∇pdxdt, stabilization_bar_LSG

using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

α = 1e6
β = 1e6
ρA = 1.0
EA = 1.0
a = 2.0
b = sqrt(2.0)
l = 4.0
c = (EA/ρA)^0.5
φ(x,y) = a*b*sin(π*x/l)*sin(π*y/l)*π
𝑢(x,y,t) = a*sin.(π.*x/l)*sin.(π.*y/l)*sin.(b*π.*t)
# ∂u∂t(x,t) = (-π.*a.*c/l)*sin.(π.*a.*c*t/l).*sin.(π.*x/l)
# ∂u∂x(x,t) = (π./l)*cos.(π.*a.*c*t/l).*cos.(π.*x/l)

const to = TimerOutput()
gmsh.initialize()
@timeit to "open msh file" gmsh.open("./msh/3D_square/4.msh")

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


@timeit to "calculate ∫∫∇q∇pdxdt" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :EA=>(x,y,z)->EA, :ρA=>(x,y,z)->ρA)
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
    @timeit to "get elements" elements_5 = getElements(nodes, entities["Γ⁵"])
    @timeit to "get elements" elements_6 = getElements(nodes, entities["Γ⁶"])
    prescribe!(elements_5,:t=>(x,y,z)->φ(x))
    prescribe!(elements_1,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_2,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_3,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_4,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_5,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    prescribe!(elements_6,:g=>(x,y,z)->0.0, :α=>(x,y,z)->α)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    @timeit to "calculate shape functions" set𝝭!(elements_5)
    @timeit to "calculate shape functions" set𝝭!(elements_6)
    𝑓 = ∫vtdΓ=>elements_5
    𝑎ᵅ = ∫vgdΓ=>elements_1∪elements_2∪elements_3∪elements_4∪elements_5
    𝑎ᵝ = ∫vgdΓ=>elements_1∪elements_2∪elements_3∪elements_4∪elements_6
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

# gmsh.finalize()

fig = Figure(resolution=(1200, 800))
ax = Axis3(fig[1,1], title="3D 位移可视化", xlabel="X", ylabel="Y", zlabel="Z")

xs = [node.x for node in nodes]
ys = [node.y for node in nodes]
zs = [node.z for node in nodes]
# ds = [node.d for node in nodes]
us = 𝑢(xs,ys,zs)

meshscatter!(
    ax, xs, ys, zs; 
    color = us,
    markersize = 0.06,
    colormap = :viridis,
    colorrange = extrema(ds)
)

Colorbar(fig[1,2], limits=extrema(us), colormap=:viridis, label="位移值")
fig

# fig = Figure()
# ax1 = Axis3(fig[1,1])
# # ax2 = Axis3(fig[1,2])

# xs = zeros(nₚ)
# ys = zeros(nₚ)
# zs = zeros(nₚ)
# ds = zeros(nₚ)
# # δds = zeros(nₚ)
# # es = zeros(nₚ)
# us = zeros(nₚ)
# # qs = zeros(nₚ)
# # as = zeros(nₚ)
# for (i, node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     z = node.z
#     us[i] = 𝑢(x,y)
#     # qs[i] = ∂u∂t(x,y)
#     # as[i] = ∂²u∂t²(x,y)
# end
# for (i,node) in enumerate(nodes)
#     xs[i] = node.x
#     ys[i] = node.y
#     zs[i] = node.z
#     ds[i] = node.d
#     # δds[i] = node.δd
#     # es[i] = ds[i] - us[i]
# end
# face = zeros(nₑ,6)
# for (i,elm) in enumerate(elements)
#     face[i,:] .= [x.𝐼 for x in elm.𝓒]
# end

# # mesh!(ax,xs,ys,face,color=zs)
# # meshscatter!(ax1,xs,ys,es,color=es,markersize = 0.1)
# meshscatter!(ax1,xs,ys,zs,ds,color=ds,markersize = 0.06)
# # meshscatter!(ax1,xs,ys,us,color=us,markersize = 0.1)
# # meshscatter!(ax2,xs,ys,δds,color=δds,markersize = 0.06)
# fig

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

# colors = zeros(nᵤ)
# 𝗠 = zeros(10)
# for (i,node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     z = node.z
#     indices = sp(x,y,z)
#     ni = length(indices)
#     𝓒 = [nodes_p[i] for i in indices]
#     data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[z]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
#     ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
#     𝓖 = [ξ]
#     a = type(𝓒,𝓖)
#     set𝝭!(a)
#     p = 0.0
#     N = ξ[:𝝭]
#     for (k,xₖ) in enumerate(𝓒)
#         p += N[k]*xₖ.p
#     end
#     colors[i] = p
# end
# α = 1.0
# points = [[node.x+α*node.u₁ for node in nodes]';[node.y+α*node.u₂ for node in nodes]';[node.z+α*node.u₃ for node in nodes]']
# # cells = [MeshCell(VTKCellTypes.VTK_TETRA,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
# vtk_grid("./vtk/block_"*poly*"_"*string(ndiv)*"_"*string(nₚ),points,cells) do vtk
#     vtk["u"] = (𝑢₁,𝑢₂,𝑢₃)
#     vtk["𝑝"] = colors
# end

