module TestMomentumGraphs
using DirectedHalfEdgeGraphs
using MomentumGraphs
using Test
using Catlab.CategoricalAlgebra.CSets
using MomentumGraphs.Form
using Symbolics
using Catlab
using Catlab.Theories, Catlab.CategoricalAlgebra

import Symbolics: variables, value
using Catlab.Graphs.BasicGraphs
using Catlab.Graphs.GraphGenerators
using Catlab.CategoricalAlgebra

import Catlab.Theories: src, tgt

using MomentumGraphs.FieldGraphs

FieldGraph()
qd = qDiagram(
  [-2, 1, 3, -1, -4, 2, 4, -3],
  [0, 2, 4, 0, 0, 1, 3, 0],
  [1, 1, 1, 1, 2, 2, 2, 2],
  [:phi1c, :photon, :photon, :phi1, :phi2c, :photon, :photon, :phi2],
  ID=1,
  pre_factor=(+1) * 1 // 2,
  nprops=2,
  nloops=1,
  nin=2,
  nout=2
)

mg = MomentumGraph([-2, 1, 3, -1, -4, 2, 4, -3],
  [0, 2, 4, 0, 0, 1, 3, 0],
  [1, 1, 1, 1, 2, 2, 2, 2])
momentum_equations_solved(mg)
dangling_edges(mg)
vcat((half_edge_pairs(mg)[1:2])...)
repeat(FVector.(value.(variables(:q, 1:5))), inner=2)

fg = FieldGraph([-2, 1, 3, -1, -4, 2, 4, -3],
  [0, 2, 4, 0, 0, 1, 3, 0],
  [1, 1, 1, 1, 2, 2, 2, 2],
  [:phi1c, :photon, :photon, :phi1, :phi2c, :photon, :photon, :phi2], dualDict=Dict(:phi1 => :phi1c, :phi2 => :phi2c, :photon => :photon), massDict=Dict(:photon => 0, :phi1 => 1, :phi2 => 2, :phi1c => 1, :phi2c => 2))
to_graphviz(fg,field_colors=Dict(:photon => "black", :phi1 => "blue", :phi2 => "red", :phi1c => "blue", :phi2c => "red"))


Hto, Hfrom, ins = to_cat_format([-2, 1, 3, -1, -4, 2, 4, -3],
  [0, 2, 4, 0, 0, 1, 3, 0])
_f(a::Real, b::Int=1) = a + 2b
f(a::Int, b::Int=0, c=4) = c - _f(a, b)

nickel_index(call_nauty(qd.g).cset)
to_graphviz(qd.g)
include("julTest2l.jl")
using StatsBase

import Catlab.Graphics.Graphviz: Graph,run_graphviz

function Base.show(io::IO, ::MIME"image/svg+xml", graph::Graph)
  run_graphviz(io, graph, format="svgz")
end

import Catlab.Theories: ⊕
function ⊕(a, b...)
  return ⊕(a, ⊕(b...))
end
⊕([qd.g, qd.g, qd.g]...)
unique!(d -> d.nickel_index, diags)
diags
(d -> d.nickel_index).(diags)
#(sample(diags,    253 , replace=false))...
alldiags=(⊕(((x -> x.g).(diags))...))
to_graphviz(alldiags,field_colors=Dict(:photon => "black", :phi1 => "blue", :phi2 => "red", :phi1c => "deepskyblue", :phi2c => "indianred"))
(x -> x.g).(diags)


FVector()
p = FVector(:p)
p(FIndex(:μ))
FVector.(value.(variables(:p, 1:4)))
FSymbol(:aμ)


mmg = FieldGraph(5)

@test add_vertex!(mmg) == 6

#add_edges!(mmg, [1,2,4,3], [2,4,3,1])
@test add_edges!(mmg, [1, 1, 1, 2, 2, 4, 3, 5, 6], [2, 5, 3, 5, 4, 6, 4, 6, 5]) == 1:18
@test add_dangling_edges!(mmg, [1, 2, 3, 4], dirs=[true, true, false, false]) == 19:22

mmg
mass(mmg, 5)
join(string.(momentum(mmg, (dangling_edges(mmg)))), " ")
zip((half_edge_pairs(mmg)))
to_graphviz(mmg)
spanningTree = subtree(mmg, dfs_parents(mmg, 4, all_neighbors, by=(x -> Base.sort(x, by=e -> mass(mmg, e), rev=true))))
to_graphviz(spanningTree, node_labels=true, edge_labels=true)

smmg = Subobject(mmg, H=half_edges(mmg), V=vertices(mmg))
to_graphviz(smmg \ spanningTree)
smmg \ spanningTree |> hom |> dom


momentum_equations(mmg)
momentum_equations_solved(mmg)


using Catlab.Graphs.GraphAlgorithms
connected_component_projection_bfs(mmg)
connected_component_projection(mmg)


mg = MomentumGraph(5)
add_edges!(mg, [1, 2, 1, 4], [2, 3, 3, 5])
connected_components(mg)
(connected_component_projection_bfs(mg))
nv(mg)


g = DirectedHalfEdgeGraph(5)

half_edges(g)
half_edges(g::AbstractDirectedHalfEdgeGraph; by::Function=(x -> sort(g, x))) = by(parts(g, :H))

half_edges(mg)
nickel_index(mmg)
(isomorphisms(mmg, mmg, iso=true))
Base.sort(all_neighbors(mmg, 3))
to_graphviz(mmg)
half_edges(mmg, 3)
mg
in_edges(mg)
dangling_edges(mg)
momentum_equations(mg)
@testset "MomentumGraphs.jl" begin
  # Write your tests here.
end
using CSetAutomorphisms

nickel_index(call_nauty(mmg).cset)



end