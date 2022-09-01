

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

function qDiagram(H::AbstractVector{Int},H′::AbstractVector{Int},vs::AbstractVector{Int},fields::AbstractVector{Symbol};ID::Int,pre_factor::Rational,nprops::Int,nloops::Int,nin::Int,nout::Int)
  g=FieldGraph(H,H′,vs,fields,dualDict=Dict(:phi1=>:phi1c,:phi2=>:phi2c,:photon=>:photon), massDict=Dict(:photon=>0,:phi1=>1,:phi2=>2,:phi1c=>1,:phi2c=>2))
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

struct Propagator <: IsProp 
  field::Symbol
  rule::Function
end

struct Theory 
  vertexrules::AbstractVector{Vertex}
  propagatorrules::AbstractVector{Propagator}
end

