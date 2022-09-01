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
using CSetAutomorphisms
using MomentumGraphs.FieldGraphs

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


to_graphviz(qd.g,field_colors=Dict(:photon => "black", :phi1 => "blue", :phi2 => "red", :phi1c => "blue", :phi2c => "red"))

nickel_index(call_nauty(qd.g).cset)

include("QGRAFjl/julTest2l.jl")
using StatsBase

import Catlab.Graphics.Graphviz: Graph,run_graphviz

function Base.show(io::IO, ::MIME"image/svg+xml", graph::Graph)
  run_graphviz(io, graph, format="svg")
end

import Catlab.Theories: ⊕
function ⊕(a, b...)
  return ⊕(a, ⊕(b...))
end

unique!(d -> d.nickel_index, diags)
diags
(d -> d.nickel_index).(diags)
#(sample(diags,    253 , replace=false))...
alldiags=(⊕(((x -> x.g).(diags))...))
to_graphviz(alldiags,field_colors=Dict(:photon => "black", :phi1 => "blue", :phi2 => "red", :phi1c => "deepskyblue", :phi2c => "indianred"))
(x -> x.g).(diags)





end