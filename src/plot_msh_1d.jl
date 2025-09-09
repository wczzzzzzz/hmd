using ApproxOperator, GLMakie

import Gmsh: gmsh

gmsh.initialize()
# gmsh.open("./msh/bar/bar_20.msh")
gmsh.open("./msh/bar/bar_un_20.msh")

entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes, entities["Ω"])
elements["Γ¹"] = getElements(nodes, entities["Γ¹"])
elements["Γ²"] = getElements(nodes, entities["Γ²"])
elements["∂Ω"] = elements["Γ¹"] ∪ elements["Γ²"]

f = Figure()

ax = Axis(f[1, 1], 
    aspect = DataAspect(), 
    xlabel = "",           
    ylabel = "",           
    xticksvisible = false,  
    yticksvisible = false, 
    xticklabelsvisible = false, 
    yticklabelsvisible = false 
)
hidespines!(ax) 

x_nodes = [node.x for node in nodes]
y_nodes = [node.y for node in nodes]
scatter!(ax, x_nodes, y_nodes,
    marker = :circle,
    markersize = 5.0,
    color = :black
)

for elm in elements["Ω"]
    x = [node.x for node in elm.𝓒]
    y = [node.y for node in elm.𝓒]
    lines!(ax, x, y,
        linewidth = 1.0,
        color = :blue
    )
end

for elm in elements["∂Ω"]
    x = [node.x for node in elm.𝓒]
    y = [node.y for node in elm.𝓒]
    scatter!(ax, x, y,
        marker = :star,
        markersize = 4.0,
        color = :red
    )
end

display(f)
# save("./fig/非均布弹簧小车/prog=1.05网格图/n=320.png",f)
