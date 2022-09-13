module Diagrams

abstract type AbstractDiagram end

struct qDiagram{G}
  ID::Int
  pre_factor::Rational
  nloops::Int
  nprops::Int
  nin::Int
  nout::Int
  g::G
  nickel_index::String 
end
using CSetAutomorphisms
using .FieldGraphs

function qDiagram(;
  H::AbstractVector{Int},
  Hdual::AbstractVector{Int},
  vertex::AbstractVector{Int},
  field::AbstractVector{Symbol},
  ID::Int,
  pre_factor::Rational,
  nprops::Int,
  nloops::Int,
  nin::Int,
  nout::Int,
  dualDict,
  massDict)
  g=FieldGraph(H,Hdual,vertex,field;dualDict,massDict)
  return qDiagram{FieldGraph}(ID,pre_factor,nprops,nloops,nin,nout,g,nickel_index(call_nauty(g).cset))
end

abstract type FeynmanRule end
abstract type IsVertex <:FeynmanRule end
abstract type IsProp <:FeynmanRule end

struct Vertex <: IsVertex 
  degree::Int
  fields::AbstractSet{Symbol}
  rule::Function  
end

seagull_rule = function(g, v)
  
  return g.H[v].degree == 1
end

struct Propagator <: IsProp 
  field::Symbol
  rule::Function
end

struct Theory 
  vertexrules::AbstractVector{Vertex}
  propagatorrules::AbstractVector{Propagator}
end

end