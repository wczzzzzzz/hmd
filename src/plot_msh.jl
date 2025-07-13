
using ApproxOperator, GLMakie

import Gmsh: gmsh

ndiv = 5
gmsh.initialize()
gmsh.open("./msh/Non-uniform/619/Tri6_2.0_"*string(ndiv)*".msh")
# gmsh.open("./msh/Non-uniform/RefineMesh_1.0/Tri6_"*string(ndiv)*".msh")
# gmsh.open("./msh/Non-uniform/局部加密/Tri3_"*string(ndiv)*".msh")
# gmsh.open("./msh/b=2/Tri3反向"*string(ndiv)*".msh")
# gmsh.open("./msh/square/square_"*string(ndiv)*".msh")
# gmsh.open("./msh/b=2/Tri3非均布"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate_"*string(ndiv)*".msh")

entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes,entities["Ω"])
elements["Γ¹"] = getElements(nodes,entities["Γ¹"])
elements["Γ²"] = getElements(nodes,entities["Γ²"])
elements["Γ³"] = getElements(nodes,entities["Γ³"])
elements["Γ⁴"] = getElements(nodes,entities["Γ⁴"])
# elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"])
# elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"])
# elements["∂Ω"] = elements["Γᵍ"]∪elements["Γᵗ"]
# elements["Γᵉ"] = getElements(nodes,entities["Γᵉ"])
elements["∂Ω"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]

# gmsh.finalize()

f = Figure()

# axis
ax = Axis3(f[1, 1], perspectiveness = 0.0, aspect = :data, azimuth = -0.5*pi, elevation = 0.5*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.))
hidespines!(ax)
hidedecorations!(ax)

x =  nodes.x
y = nodes.y
z = 0
ps = Point3f.(x,y,z)
scatter!(ps, 
    marker=:circle,
    markersize = 8.0,
    color = :black
)

elements
for elm in elements["Ω"]
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]
    # x = [x.x for x in elm.𝓒[[1,2,3,4]]]
    # y = [x.y for x in elm.𝓒[[1,2,3,4]]]

    lines!(x,y, linewidth = 1.5, color = :black)
end

# # boundaries
for elm in elements["∂Ω"]
    ξ¹ = [x.x for x in elm.𝓒]
    ξ² = [x.y for x in elm.𝓒]
    x =  [x.x for x in elm.𝓒]
    y =  [x.y for x in elm.𝓒]
    lines!(x,y,linewidth = 1.5, color = :black)
end

# save("./fig/square_"*string(ndiv)*".png",f)
# save("./fig/三角形节点网格/Tri3非均布_Rf_1.0_5"*string(ndiv)*".png",f, px_per_unit = 2.0)
# save("./fig/三角形节点网格/Tri3_square_"*string(ndiv)*".png",f)
# save("./fig/三角形节点网格/Tri6非均布_Rf_1.0_3"*string(ndiv)*".png",f, px_per_unit = 2.0)

# save("./fig/619测试/Tri6_网格图非均布2.0_"*string(ndiv)*".png",f)



f