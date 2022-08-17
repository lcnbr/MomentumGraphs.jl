module MomentumGraphs

module Form
using SymbolicUtils
import SymbolicUtils: to_symbolic
import Symbolics: variable,value,getname
export FIndex,FSymbol,FVector
import Base: +, -, *


struct FSymbol
  symbol::SymbolicUtils.Sym
end

function ⋅(a, b) end
@register ⋅(a, b) 

SymbolicUtils.to_symbolic(s::FSymbol) = s.symbol

⋅(a::FSymbol, b::FSymbol) = FSymbol(a.symbol ⋅ b.symbol)
+(a::FSymbol , b::FSymbol ) = FSymbol(a.symbol + b.symbol)
*(a::FSymbol , b::FSymbol ) = FSymbol(a.symbol * b.symbol)

FSymbol(s::Symbol)=FSymbol(value(variable(s)))

Base.show(io::IO, s::FSymbol) = print(io,  s.symbol)
Base.show(io::IO, ::MIME"text/plain", s::FSymbol) =
           print(io, "FORM Symbol: ", s)
Base.show(io::IO, ::MIME"text/FORM", s::FSymbol) =
           print(io, "S", s,";")


struct FIndex 
  dimension::Union{Int,FSymbol}
  ID::Symbol
  generated::Bool
  function  FIndex(dimension,ID)
    new(dimension,ID,is_gensym(ID))
  end
end

function FIndex(ID=gensym(:dummy);dimension=4,)
  FIndex(dimension,ID)
end

function Base.show(io::IO, i::FIndex)
  if(i.generated)
    print(io, "μ",  subscript(gensym_to_num(i.ID)))
  else
    print(io, i.ID)
  end
end

Base.show(io::IO, ::MIME"text/plain", i::FIndex) =
           print(io, "FORM Index: ", i, "with dimension ", i.dimension)
Base.show(io::IO, ::MIME"text/FORM", i::FIndex) =
           print(io, "I", i,";")


mutable struct FVector 
  symbol::FSymbol
  index::FIndex
end


⋅(a::FVector , b::FVector ) = FVector(a.symbol ⋅ b.symbol,a.index)
+(a::FVector  , b::FVector ) = FVector(a.symbol + b.symbol)

function *(a::FVector  , b::FVector) 
  if a.index.ID == b.index.ID
    ⋅(a, b)
  else
    FVector(a.symbol * b.symbol, a.index)#This is clearly very wrong! 
  end
  FVector(FSymbol(symbol),FIndex(dimension=dimension))
end
FVector(a.symbol * b.symbol)

FVector(symbol;dimension=4)=FVector(FSymbol(symbol),FIndex(dimension=dimension))

FVector()=FVector(:p)

function (v::FVector)(i::FIndex)
  v.index=i
  v
end
  


Base.show(io::IO, v::FVector) = print(io, "$(v.symbol)($(v.index))",)
Base.show(io::IO, ::MIME"text/plain", v::FVector) =
           print(io, "FORM Vector: ", v)
Base.show(io::IO, ::MIME"text/FORM", v::FVector) =
           print(io, "I", v,";")

#= mutable struct FTensor 
  symbol::FSymbol
  indices::Vector{FIndex}
end
          
          
⋅(a::FTensor , b::FTensor ) = FTensor(a.symbol ⋅ b.symbol,a.index)
+(a::FTensor  , b::FTensor ) = FTensor(a.symbol + b.symbol)
          
function *(a::FTensor  , b::FTensor) 
  if a.index.ID == b.index.ID
    ⋅(a, b)
  else
    FTensor(a.symbol ⋅ b.symbol, a.index)
  end
  FTensor(FSymbol(symbol),FIndex(dimension=dimension))
end
FTensor(a.symbol * b.symbol)

FTensor(symbol;dimension=4)=FTensor(FSymbol(symbol),FIndex(dimension=dimension))

FTensor()=FTensor(:p)

function (v::FTensor)(i::FIndex)
  v.index=i
  v
end =#


function gennum(a)
  gensym_to_num(gensym(a))
end

function gensym_to_num(symbol)
  @assert is_gensym(symbol) "Symbol has not been generated"
  parse(Int,match(r".*#\K\d*",string( symbol)).match)-291
end

function is_gensym(symbol::SymbolicUtils.Sym)
  is_gensym(getname(symbol))
end

function is_gensym(symbol::Symbol)
  occursin(r"^##",string( symbol))
end

subscript(i) = join(Char(0x2080 + d) for d in reverse!(digits(i)))

end

using .Form
import Symbolics: variables,value
using DirectedHalfEdgeGraphs


export SchMomentumGraph,AbstractMomentumGraph,MomentumGraph,sink,source,in_edges,out_edges,half_edge_pairs,add_half_edge_pairs!,add_dangling_edge!,add_dangling_edges!,to_graphviz,to_graphviz_property_graph,momentum

using Catlab
using Base: @invoke
using Catlab.CategoricalAlgebra.CSets,Catlab.Graphics,Catlab.Graphs,Catlab.Graphics.GraphvizGraphs

import Catlab.Graphs.BasicGraphs: add_dangling_edges!, add_dangling_edge!, half_edge_pairs,add_half_edge_pairs!
import Catlab.Graphics.GraphvizGraphs: to_graphviz, to_graphviz_property_graph

@present SchMomentumGraph <: SchDirectedHalfEdgeGraph begin
  name::AttrType
  momentum::Attr(H, name)
  compose(inv,momentum)==momentum
end

@abstract_acset_type AbstractMomentumGraph <: AbstractDirectedHalfEdgeGraph


@acset_type MomentumGraphGeneric(SchMomentumGraph, index=[:inv, :vertex, :sink]) <: AbstractMomentumGraph

MomentumGraph=MomentumGraphGeneric{Bool,FVector}

momentum(g::AbstractMomentumGraph, args...) = subpart(g, args..., :momentum)


function add_half_edge_pairs!(g::AbstractMomentumGraph, srcs::AbstractVector{Int},
  tgts::AbstractVector{Int}; kw...)

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
  first(outs):last(ins)
end

function add_dangling_edge!(g::AbstractMomentumGraph, v::Int; dir=true, kw...)
  neIn=length(dangling_edges(g))
  H=add_part!(g, :H; vertex=v, inv=nparts(g,:H)+1,sink=dir,momentum=FVector(value(variables(:p,neIn+1)[1])),kw...)
end

function add_dangling_edges!(g::AbstractMomentumGraph, vs::AbstractVector{Int}; dirs::AbstractVector{Bool},kw...)
  neIn=length(dangling_edges(g))
  n, k = length(vs), nparts(g, :H)
  H=add_parts!(g, :H, n; vertex=vs, inv=(k+1):(k+n),sink=dirs, momentum=FVector.(value.(variables(:p,(neIn+1):neIn+n))),kw...)
end

function momentum_equations(g::AbstractMomentumGraph)
  equations(g, :momentum)
end
  
end

end






