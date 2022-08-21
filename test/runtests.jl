module TestMomentumGraphs
using DirectedHalfEdgeGraphs
using MomentumGraphs
using Test
using Catlab.CategoricalAlgebra.CSets
using MomentumGraphs.Form

using Catlab

import Symbolics: variables,value
using Catlab.Graphs.BasicGraphs
using Catlab.Graphs.GraphGenerators

import Catlab.Theories: src, tgt

MassiveMomentumGraph()

MomentumGraphs.MassiveMomentumGraph(3)

FVector()
p=FVector(:p)
p(FIndex(:μ))
FVector.(value.(variables(:p,1:4)))
FSymbol(:aμ)
gh = @acset DirectedHalfEdgeGraph begin
  V = 4
  H = 7
  vertex=[1,2,3,4,1,2,3]
  inv=[2,1,3,5,4,7,6]
  sink=[1,1,0,1,0,1,0].>0
end

gh=HalfEdgeGraph()
a=fVector()
b=fVector(:a,4)

add_vertices!(gh, 4)
add_edges!(gh, [1,2,3], [2,3,1])
add_dangling_edges!(gh, [5,6,3])
half_edges(gh)
subpart(gh,:inv)
gh
incident(gh,2,:vertex)
to_graphviz(gh)
dom(hom(Subobject(gh,H=1:5,V=1:4)))




g = erdos_renyi(Graph, 20, 50)
to_graphviz(g)
src(g)
mg=MomentumGraph()
add_vertices!(mg, 6)
add_edges!(mg, [1,1,1,2,2,4,4,5,6],[2,5,3,5,4,6,3,6,5] )
add_dangling_edges!(mg, [1,2,3,4],dirs=[true,false,true,false])

to_graphviz(mg)
spanningTree=subtree(mg,dfs_parents(mg,1,all_neighbors))
to_graphviz(spanningTree)
mgt=copy(mg)
rem_edge!(mg,5,6)
to_graphviz(mgt)
mg
half_edges(mg)
out_edges(mg,3)
dangling_edges(mg)
edges(mg)
edges(mg,5,6) 


to_graphviz(mg)
mg

connected_component_projection_bfs(mg)
connected_component_projection(mg)
connected_components(mg)
nv(mg)


momentum_equations(mg)
@testset "MomentumGraphs.jl" begin
    # Write your tests here.
end



end