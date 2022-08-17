module TestMomentumGraphs
using DirectedHalfEdgeGraphs
using MomentumGraphs
using Test
using Catlab.CategoricalAlgebra.CSets
using MomentumGraphs.Form
using Catlab

import Symbolics: variables,value
using Catlab.Graphs.BasicGraphs

FVector()
p=FVector(:p)
p(FIndex(:μ))
FVector.(value.(variables(:p,1:4)))
FSymbol(:aμ)
gh = @acset DirectedHalfEdgeGraph begin
  V = 4
  H = 7
  vertex= [1,2,3,4,5,6,7]
  inv=    [2,1,3,5,4,7,6]
  sink=   [1,0,0,0,1,1,0].>0

end

gh=DirectedHalfEdgeGraph()
a=fVector()
b=fVector(:a,4)

add_vertices!(gh, 4)
add_edges!(gh, [1,2,3], [2,3,1])
add_dangling_edges!(gh, [5,6,3],dirs=[true,false,true])
half_edges(gh)
subpart(gh,:inv)
incident(gh,2,:vertex)
to_graphviz(gh)


mg=MomentumGraph()
add_vertices!(mg, 6)
add_edges!(mg, [1,2,3], [2,3,1])
add_dangling_edges!(mg, [5,6,3],dirs=[true,false,true])
mg
to_graphviz(mg)
half_edges(mg,4)
out_edges(mg,3)
dangling_edges(mg)
DirectedHalfEdgeGraphs.dangling_edges(mg,3)
for v ∈ vertices(mg)
  in_edges(mg,v)
  
end

(MomentumGraphs.momentum.(Ref(mg),in_edges(mg,3)))

@testset "MomentumGraphs.jl" begin
    # Write your tests here.
end
end