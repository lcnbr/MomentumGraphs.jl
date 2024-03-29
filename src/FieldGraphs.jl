module FieldGraphs

export FieldGraph, field
using DirectedHalfEdgeGraphs
using ...Form
using ...Fields
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.CSets, Catlab.Graphics, Catlab.Graphs, Catlab.Graphics.GraphvizGraphs


import DirectedHalfEdgeGraphs: add_dangling_edges!, add_dangling_edge!,
  half_edge_pairs, add_half_edge_pairs!,
  add_half_edges!, nickel_index, sort, to_graphviz, to_graphviz_property_graph
import ...MassiveMomentumGraphs: SchMassiveMomentumGraph, AbstractMassiveMomentumGraph, momentum, mass

import Symbolics: variables, value, solve_for

@present SchFieldGraph <: SchMassiveMomentumGraph begin
  Field::AttrType
  field::Attr(H, Field)
end


@abstract_acset_type AbstractFieldGraph <: AbstractMassiveMomentumGraph


@acset_type FieldGraphGeneric(SchFieldGraph, index=[:inv, :vertex, :sink, :momentum, :indep, :mass, :field]) <: AbstractFieldGraph

FieldGraph = FieldGraphGeneric{Bool,FVector,Real,AbstractField}


mass(f::AbstractField) = Fields.mass(f)
mass(::Type{T}) where {T<:AbstractField} = Fields.mass(T)


import Fields: field
field(g::AbstractFieldGraph, args...) = subpart(g, args..., :field)
sort(g::AbstractFieldGraph, args...) = sort(args..., by=e -> mass(g, e))

function add_half_edge_pairs!(g::AbstractFieldGraph, srcs,
  tgts;
  indeps=falses(length(srcs)), fieldtypes=[ScalarField{:s} for i in 1:length(srcs)], kw...)

  @assert (n = length(srcs)) == length(tgts) == length(indeps) == length(fieldtypes)

  neIn = length(first(half_edge_pairs(g)))

  momenta = FVector.(value.(variables(:q, (neIn+1):neIn+n)))

  outs = add_parts!(g, :H, n; vertex=srcs, kw...)
  ins = add_parts!(g, :H, n; vertex=tgts, kw...)
  set_subpart!(g, outs, :inv, ins)
  set_subpart!(g, outs, :sink, falses(n))
  set_subpart!(g, ins, :inv, outs)
  set_subpart!(g, ins, :sink, trues(n))
  set_subpart!(g, outs, :momentum, momenta)
  set_subpart!(g, ins, :momentum, momenta)
  set_subpart!(g, outs, :indep, indeps)
  set_subpart!(g, ins, :indep, indeps)
  set_subpart!(g, outs, :mass, Vector{Int}(.!ismassless.(fieldtypes)))
  set_subpart!(g, ins, :mass, Vector{Int}(.!ismassless.(fieldtypes)))

  set_subpart!(g, outs, :field, map.(fieldtypes, outs))
  set_subpart!(g, ins, :field, map.(adjoint.(fieldtypes), ins))

  first(outs):last(ins)
end

function add_dangling_edge!(g::AbstractFieldGraph, v::Int;
  dir=true, indep=true,
  momentum=FVector(value(variables(:p, length(dangling_edges(g)) + 1)[1])),
  fieldtype::Type{AbstractField}=1:ScalarField{:s}, kw...)

  h = add_part!(g, :H; vertex=v, inv=nparts(g, :H) + 1, sink=dir)
  set_subparts!(g, h; momentum=momentum, indep=indep, mass=Int(!ismassless(fieldtype)))
  set_subparts!(g, h; field=dir ? fieldtype(h) : adjoint(fieldtype(h)), kw...)
  h
end

function add_dangling_edges!(g::AbstractFieldGraph, vs::AbstractVector{Int};
  dirs::AbstractVector{Bool}=trues(length(vs)), indeps=trues(length(vs)),
  momenta=FVector.(value.(variables(:p, (length(dangling_edges(g))+1):length(dangling_edges(g))+length(vs)))),
  fieldtypes::AbstractVector{DataType}=[ScalarField{:s} for i in 1:length(vs)], kw...)

  n, k = length(vs), nparts(g, :H)
  hs = add_parts!(g, :H, n; vertex=vs, inv=(k+1):(k+n), sink=dirs)
  set_subparts!(g, hs; momentum=momenta, indep=indeps, mass=Vector{Int}(.!ismassless.(fieldtypes)))
  set_subparts!(g, hs; field=[field = dirs[i] ? fieldtype(hs[i]) : adjoint(fieldtype(hs[i])) for (i, fieldtype) ∈ enumerate(fieldtypes)], kw...)
  hs
end


function add_half_edges!(g::AbstractFieldGraph, invs::Vector{Int}, vertex::Vector{Int}, sinks::AbstractVector{Bool}=falses(length(inv)),
  fieldtypes::AbstractVector{DataType}=[ScalarField{:s} for i in 1:length(invs)];
  strict=false, kw...)

  niIn = length(first(half_edge_pairs(g)))
  niExt = length(dangling_edges(g))
  he = _add_half_edges!(g, invs, vertex, sinks; strict, kw...)
  set_subpart!(g, he, :field, map.(fieldtypes, he))
  if strict
    for h in he

      @show typeof(field(g, h)) sink(g, h) adjoint(typeof(field(g, inv(g, h))))
      @assert inv(g, h) == h || typeof(field(g, h)) == adjoint(typeof(field(g, inv(g, h))))

      @assert isselfadjoint(field(g, h)) || isadjoint(field(g, h)) == sink(g, h)
    end
  else
    for h in he
      if sink(g, h)
        set_subpart!(g, h, :field, adjoint(basefield(field(g, h))))
      else
        set_subpart!(g, h, :field, basefield(field(g, h)))
      end
    end
  end


  set_subpart!(g, he, :mass, Vector{Int}(.!ismassless.(fieldtypes)))
  nfIn = length(first(half_edge_pairs(g)))
  nfExt = length(dangling_edges(g))

  innermomenta = FVector.(value.(variables(:q, (niIn+1):nfIn)))
  inh = half_edge_pairs(g)

  set_subpart!(g, inh[1][(niIn+1):nfIn], :momentum, innermomenta)
  set_subpart!(g, inh[2][(niIn+1):nfIn], :momentum, innermomenta)


  outermomenta = FVector.(value.(variables(:p, (niExt+1):nfExt)))
  exth = dangling_edges(g)[(niExt+1):nfExt]
  set_subpart!(g, exth, :momentum, outermomenta)
  he
end



function nickel_index(g::AbstractFieldGraph)
  index = edge_index(g)
  index *= " : "
  for v ∈ vertices(g)

    momenta = join(string.(momentum(g, (dangling_edges(g, v)))), " ")
    index *= string("|", momenta)
    ns = Base.sort(all_neighbors(g, v))
    newns = ns[ns.>v]
    index *= string(join(string.("M", mass(g, newns), " ")))

  end

  index
end


function default_node_attrs(labels::Union{Symbol,Bool})
  shape = labels isa Symbol ? "ellipse" : (labels ? "circle" : "point")
  Dict(:shape => shape, :width => "0.05", :height => "0.05", :margin => "0")
end
node_label(g, name::Symbol, v::Int) = Dict(:label => string(g[v, name]))
node_label(g, labels::Bool, v::Int) = Dict(:label => labels ? string(v) : "")

edge_label(g, name::Symbol, e::Int) = Dict(:label => string(g[e, name]))
edge_label(g, labels::Bool, e::Int) =
  labels ? Dict(:label => string(e)) : Dict{Symbol,String}()

to_graphviz(g::AbstractFieldGraph; kw...) =
  to_graphviz(to_graphviz_property_graph(g; kw...))

function to_graphviz_property_graph(g::AbstractFieldGraph;
  prog::AbstractString="neato", graph_attrs::AbstractDict=Dict(),
  node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(),
  node_labels::Union{Symbol,Bool}=false, edge_labels::Union{Symbol,Bool}=false)
  pg = PropertyGraph{Any}(; prog=prog,
    graph=graph_attrs,
    node=merge!(default_node_attrs(node_labels), node_attrs),
    edge=merge!(Dict(:arrowsize => "0.5"), edge_attrs)
  )
  for v in vertices(g)
    add_vertex!(pg, label=node_labels ? string(v) : "")
  end
  for h in dangling_edges(g)
    # Dangling half-edge: add an invisible vertex.
    v = add_vertex!(pg, style="invis", shape="none", label="")
    if sink(g, h)
      e′ = add_edge!(pg, v, vertex(g, h), penwidth=string(mass(g, h) + 1), color=color(field(g, h)))
    else
      e′ = add_edge!(pg, vertex(g, h), v, penwidth=string(mass(g, h) + 1), color=color(field(g, h)))
    end
    set_eprops!(pg, e′, edge_label(g, edge_labels, h))
  end
  for (source, sink) ∈ zip(half_edge_pairs(g)...)
    (src, tgt) = (vertex(g, source), vertex(g, sink))
    e = add_edge!(pg, src, tgt, penwidth=string(mass(g, source) + 1), color=string(color(field(g, source)), ";0.5:", color(field(g, sink))))
    set_eprops!(pg, e, edge_label(g, edge_labels, e))
  end
  pg
end



end