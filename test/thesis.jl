module Thesis
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
using CSetAutomorphisms
using MomentumGraphs.FieldGraphs


dualDict=Dict(:phi1=>:phi1c,:phi2=>:phi2c,:photon=>:photon)
massDict=Dict(:photon=>0,:phi1=>1,:phi2=>2,:phi1c=>1,:phi2c=>2)
field_colors=Dict(:photon => "black", :phi1 => "blue", :phi2 => "red", :phi1c => "deepskyblue", :phi2c => "indianred")
diags=include("QGRAFjl/1lSQED.jl")

qDiags=[qDiagram(;diag...,dualDict,massDict) for diag in diags]

import Catlab.Theories: ⊕
function ⊕(a, b...)
  return ⊕(a, ⊕(b...))
end

unique!(d -> d.nickel_index, qDiags)
qDiags
(d -> d.nickel_index).(qDiags)
#(sample(qDiags,    253 , replace=false))...
alldiags=(⊕(((x -> x.g).(qDiags))...))
to_graphviz(alldiags,field_colors=Dict(:photon => "black", :phi1 => "blue", :phi2 => "red", :phi1c => "deepskyblue", :phi2c => "indianred"))


qGraphs=(x -> x.g).(qDiags)
to_graphviz(qGraphs[13];field_colors)
momentum_equations_solved.(qGraphs;by=massSort)

massSort=(x -> Base.sort(x, by=e -> mass(box.g, e), rev=true))

box=qDiags[12]
to_graphviz(box.g;field_colors)
momentum_equations(box.g) 
momentum_equations_solved(box.g;by=massSort)
box.g
spanningTree = subtree(box.g, dfs_parents(box.g, 4, all_neighbors, by=massSort))
to_graphviz(spanningTree, node_labels=true, edge_labels=true; field_colors)
sbox = Subobject(box.g, H=half_edges(box.g), V=vertices(box.g))
to_graphviz(sbox \ spanningTree;field_colors)

crossBox =qDiags[13]
to_graphviz(crossBox.g;field_colors)


end