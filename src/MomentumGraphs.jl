module MomentumGraphs


include("Form.jl")

using .Form
import Symbolics: variables,value
using DirectedHalfEdgeGraphs


export SchMomentumGraph,AbstractMomentumGraph,MomentumGraph,sink,source,in_edges,out_edges,half_edge_pairs,add_half_edge_pairs!,add_dangling_edge!,add_dangling_edges!,to_graphviz,to_graphviz_property_graph,momentum,momentum_equation,momentum_equations,SchMassiveMomentumGraph,AbstractMassiveMomentumGraph,MassiveMomentumGraph,mass

using Catlab
using Base: @invoke
using Catlab.CategoricalAlgebra.CSets,Catlab.Graphics,Catlab.Graphs,Catlab.Graphics.GraphvizGraphs

import Catlab.Graphs.BasicGraphs: add_dangling_edges!, add_dangling_edge!, half_edge_pairs,add_half_edge_pairs!
import Catlab.Graphics.GraphvizGraphs: to_graphviz, to_graphviz_property_graph

@present SchMomentumGraph <: SchDirectedHalfEdgeGraph begin
  Vec::AttrType
  momentum::Attr(H, Vec)
  indep::Attr(H, Truth)
  compose(inv,momentum)==momentum
end

@abstract_acset_type AbstractMomentumGraph <: AbstractDirectedHalfEdgeGraph


@acset_type MomentumGraphGeneric(SchMomentumGraph, index=[:inv, :vertex, :sink,:momentum,:indep]) <: AbstractMomentumGraph

MomentumGraph=MomentumGraphGeneric{Bool,FVector}


#******************************************************************************
# Accessors


momentum(g::AbstractMomentumGraph, args...) = subpart(g, args..., :momentum)
function momentum_equation(g::AbstractMomentumGraph,vertex)
  ins=(p->p.symbol).((p->p.symbol).(momentum.(Ref(g),in_edges(g,vertex))))
  outs=(p->p.symbol).((p->p.symbol).(momentum.(Ref(g),out_edges(g,vertex))))
  length(ins)>0 ? ins=sum(ins) : ins=0
  length(outs)>0 ? outs=sum(outs) : outs=0
  ins~outs
end
function momentum_equations(g::AbstractMomentumGraph)

  [momentum_equation(g,v) for v in vertices(g)]
end

#*****************************************************************************
# Edit graph

function add_half_edge_pairs!(g::AbstractMomentumGraph, srcs::AbstractVector{Int},
  tgts::AbstractVector{Int}; indeps=falses(length(srcs)),kw...)

  @assert (n = length(srcs)) == length(tgts)
  
  neIn=length(first(half_edge_pairs(g)))

  momenta=FVector.(value.(variables(:q,(neIn+1):neIn+n)))

  outs = add_parts!(g, :H, n; vertex=srcs, kw...)
  ins = add_parts!(g, :H, n; vertex=tgts, kw...)
  set_subpart!(g, outs, :inv, ins)
  set_subpart!(g, outs, :sink, falses(n))
  set_subpart!(g, ins, :inv, outs)
  set_subpart!(g, ins, :sink, trues(n))
  set_subpart!(g, outs, :momentum,momenta)
  set_subpart!(g, ins, :momentum,momenta)
  set_subpart!(g, outs, :indep, indeps)
  set_subpart!(g, ins, :indep,indeps)
  first(outs):last(ins)
end

function add_dangling_edge!(g::AbstractMomentumGraph, v::Int; dir=true,momentum=FVector(value(variables(:p,length(dangling_edges(g))+1)[1])),indep=true, kw...)
  
  H=add_part!(g, :H; vertex=v, inv=nparts(g,:H)+1,sink=dir,momentum=momentum,indep=indep,kw...)
end

function add_dangling_edges!(g::AbstractMomentumGraph, vs::AbstractVector{Int}; dirs::AbstractVector{Bool}=trues(length(vs)),momenta=FVector.(value.(variables(:p,(length(dangling_edges(g))+1):length(dangling_edges(g))+length(vs)))),indeps=trues(length(vs)),kw...)
  neIn=length(dangling_edges(g))
  n, k = length(vs), nparts(g, :H)
  H=add_parts!(g, :H, n; vertex=vs, inv=(k+1):(k+n),sink=dirs, momentum=momenta,indep=indeps,kw...)
end


export SchMassiveMomentumGraph,AbstractMassiveMomentumGraph,MassiveMomentumGraph,mass

include("MassiveMomentumGraphs.jl")
  
end
  









