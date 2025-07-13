
using ApproxOperator
import ApproxOperator.Test: cc𝝭, cc∇𝝭
using TimerOutputs 
using SparseArrays
import Gmsh: gmsh

# constant
# u(x,y,z) = 1.0
# ∂u∂x(x,y,z) = 0.0
# ∂u∂y(x,y,z) = 0.0

# linear
# u(x,y,z) = 1+2x+3y
# ∂u∂x(x,y,z) = 2.0
# ∂u∂y(x,y,z) = 3.0

# quadratic
# u(x,y,z) = 1+2x+3y+4x^2+5x*y+6y^2
# ∂u∂x(x,y,z) = 2.0+8x+5y
# ∂u∂y(x,y,z) = 3.0+12y+5x

# cubic
# u(x,y,z) = 1+2x+3y+4x^2+5x*y+6y^2+7x^3+8x^2*y+9x*y^2+10y^3
# ∂u∂x(x,y,z) = 2.0+8x+5y+21x^2+16x*y+9y^2
# ∂u∂y(x,y,z) = 3.0+12y+5x+8x^2+18x*y+30y^2

@timeit to "open msh file" gmsh.open("./msh/square/square_16.msh")
# @timeit to "open msh file" gmsh.open("./msh/Non-uniform/Tri3/4.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

@timeit to "get elements" elements = getElements(nodes, entities["Ω"])
@timeit to "get elements" elements_t, nodes_t, edges = Tri3toTriHermite(elements, nodes)
prescribe!(elements_t, :u=>u, :∂u∂x=>∂u∂x, :∂u∂y=>∂u∂y)
@timeit to "calculate shape functions" set𝝭!(elements_t)
@timeit to "calculate shape functions" set∇𝝭!(elements_t)

@timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"])
@timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"])
@timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"])
@timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"])
@timeit to "get elements" elements_1t = Seg2toSegHermite(elements_1, nodes_t, edges)
@timeit to "get elements" elements_2t = Seg2toSegHermite(elements_2, nodes_t, edges)
@timeit to "get elements" elements_3t = Seg2toSegHermite(elements_3, nodes_t, edges)
@timeit to "get elements" elements_4t = Seg2toSegHermite(elements_4, nodes_t, edges)
prescribe!(elements_1t,:u=>u)
prescribe!(elements_2t,:u=>u)
prescribe!(elements_3t,:u=>u)
prescribe!(elements_4t,:u=>u)
@timeit to "calculate shape functions" set𝝭!(elements_1t)
@timeit to "calculate shape functions" set𝝭!(elements_2t)
@timeit to "calculate shape functions" set𝝭!(elements_3t)
@timeit to "calculate shape functions" set𝝭!(elements_4t)

nₚ = length(nodes)
nₑ = length(elements_t)
nₗ = length(nodes_t) - nₚ - nₑ

# Consistency condition
d = zeros(nₚ+nₗ+nₑ)
for i in 1:nₚ
    I = nodes_t[i].𝐼
    x = nodes_t[i].x
    y = nodes_t[i].y
    z = nodes_t[i].z
    d[I] = u(x,y,z)
end
for i in 1:nₗ
    I = nodes_t[nₚ+i].𝐼
    x = nodes_t[nₚ+i].x
    y = nodes_t[nₚ+i].y
    z = nodes_t[nₚ+i].z
    s₁ = nodes_t[nₚ+i].s₁
    s₂ = nodes_t[nₚ+i].s₂
    d[nₚ+i] = ∂u∂x(x,y,z)*s₁ + ∂u∂y(x,y,z)*s₂
end
for i in 1:nₑ
    I = nodes_t[nₚ+nₗ+i].𝐼
    x = nodes_t[nₚ+nₗ+i].x
    y = nodes_t[nₚ+nₗ+i].y
    z = nodes_t[nₚ+nₗ+i].z
    d[nₚ+nₗ+i] = u(x,y,z)
end
push!(nodes_t,:d=>d)

# err = cc𝝭(elements["Ωᵗ"])
# err = cc𝝭(elements["Γ₁ᵗ"])
err = cc∇𝝭(elements_t)
# err = cc𝝭(elements["Ωᵗ"][1:1])