module TestFieldGraphs
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

dualDict=Dict(:phi1 => :phi1c, :phi2 => :phi2c, :photon => :photon)
massDict=Dict(:photon => 0, :phi1 => 1, :phi2 => 2, :phi1c => 1, :phi2c => 2)
field_colors=Dict(:photon => "black", :phi1 => "blue", :phi2 => "red", :phi1c => "blue", :phi2c => "red")

fg = FieldGraph([-2, 1, 3, -1, -4, 2, 4, -3],
  [0, 2, 4, 0, 0, 1, 3, 0],
  [1, 1, 1, 1, 2, 2, 2, 2],
  [:phi1c, :photon, :photon, :phi1, :phi2c, :photon, :photon, :phi2]; dualDict, massDict)

to_graphviz(fg;field_colors)


fg = FieldGraph(5)

@test add_vertex!(fg) == 6
#add_edges!(fg, [1,2,4,3], [2,4,3,1])
@test add_edges!(fg, [1, 1, 1, 2, 2, 4, 3, 5, 6], [2, 5, 3, 5, 4, 6, 4, 6, 5]) == 1:18
@test add_dangling_edges!(fg, [1, 2, 3, 4], dirs=[true, true, false, false]) == 19:22

fg
to_graphviz(fg)

fg = FieldGraph(5)

@test add_vertex!(fg) == 6
fields=[:ϕ,:ψ,:γ,:ϕ,:ψ,:γ,:ϕ,:ψ,:γ]
dualDict=Dict(:ϕ => :ϕ̄ , :ψ => :ψ̄, :γ  => :γ)
massDict=Dict(:γ => 0, :ϕ => 1, :ψ => 2, :ϕ̄ => 1, :ψ̄ => 2)
field_colors=Dict(:γ => "black", :ψ => "blue", :ϕ => "red", :ϕ̄ => "red", :ψ̄ => "blue")

#add_edges!(fg, [1,2,4,3], [2,4,3,1])

get.(Ref(massDict), fields, missing)
get.(Ref(dualDict), fields, missing)
@test add_edges!(fg, [1, 1, 1, 2, 2, 4, 3, 5, 6], [2, 5, 3, 5, 4, 6, 4, 6, 5];fields,massDict,dualDict) == 1:18
@test add_dangling_edges!(fg, [1, 2, 3, 4]; dirs=[true, true, false, false],fields=[:ϕ,:ϕ,:ψ,:ψ],dualDict,massDict) == 19:22

fg
to_graphviz(fg;field_colors)


spanningTree = subtree(fg, dfs_parents(fg, 4, all_neighbors, by=(x -> Base.sort(x, by=e -> mass(fg, e), rev=false))))
to_graphviz(spanningTree, node_labels=true, edge_labels=true; field_colors)

sfg = Subobject(fg, H=half_edges(fg), V=vertices(fg))
to_graphviz(sfg \ spanningTree;field_colors)
sfg \ spanningTree |> hom |> dom


momentum_equations(fg)
momentum_equations_solved(fg;by=(x -> Base.sort(x, by=e -> mass(fg, e), rev=false)))


using Catlab.Graphs.GraphAlgorithms
connected_component_projection_bfs(fg)
connected_component_projection(fg)

using CSetAutomorphisms

nickel_index(call_nauty(fg).cset)



end