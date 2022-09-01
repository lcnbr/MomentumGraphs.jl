


@present SchMassiveMomentumGraph <: SchMomentumGraph begin
  Weight::AttrType
  mass::Attr(H, Weight)
end

@abstract_acset_type AbstractMassiveMomentumGraph <: AbstractMomentumGraph


@acset_type MassiveMomentumGraphGeneric(SchMassiveMomentumGraph, index=[:inv, :vertex, :sink, :momentum, :indep, :mass]) <: AbstractMassiveMomentumGraph

MassiveMomentumGraph = MassiveMomentumGraphGeneric{Bool,FVector,Real}

mass(g::AbstractMassiveMomentumGraph, args...) = subpart(g, args..., :mass)
sort(g::AbstractMassiveMomentumGraph, args...) = sort(args..., by=e -> mass(g, e))

function add_half_edge_pairs!(g::AbstractMassiveMomentumGraph, srcs::AbstractVector{Int},
  tgts::AbstractVector{Int}; indeps=falses(length(srcs)), mass=zeros(length(srcs)), kw...)

  @assert (n = length(srcs)) == length(tgts)

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
  set_subpart!(g, outs, :mass, mass)
  set_subpart!(g, ins, :mass, mass)

  first(outs):last(ins)
end

function add_dangling_edge!(g::AbstractMassiveMomentumGraph, v::Int; dir=true, momentum=FVector(value(variables(:p, length(dangling_edges(g)) + 1)[1])), indep=true, mass=0.0, kw...)

  H = add_part!(g, :H; vertex=v, inv=nparts(g, :H) + 1, sink=dir, momentum=momentum, indep=indep, mass=mass, kw...)
end

function add_dangling_edges!(g::AbstractMassiveMomentumGraph, vs::AbstractVector{Int}; dirs::AbstractVector{Bool}=trues(length(vs)), momenta=FVector.(value.(variables(:p, (length(dangling_edges(g))+1):length(dangling_edges(g))+length(vs)))), indeps=trues(length(vs)), masses=zeros(length(vs)), kw...)
  neIn = length(dangling_edges(g))
  n, k = length(vs), nparts(g, :H)
  H = add_parts!(g, :H, n; vertex=vs, inv=(k+1):(k+n), sink=dirs, momentum=momenta, indep=indeps, mass=masses, kw...)
end

#******************************************************************************
# Displaying  


function default_node_attrs(labels::Union{Symbol,Bool})
  shape = labels isa Symbol ? "ellipse" : (labels ? "circle" : "point")
  Dict(:shape => shape, :width => "0.05", :height => "0.05", :margin => "0")
end
node_label(g, name::Symbol, v::Int) = Dict(:label => string(g[v, name]))
node_label(g, labels::Bool, v::Int) = Dict(:label => labels ? string(v) : "")

edge_label(g, name::Symbol, e::Int) = Dict(:label => string(g[e, name]))
edge_label(g, labels::Bool, e::Int) =
  labels ? Dict(:label => string(e)) : Dict{Symbol,String}()

to_graphviz(g::AbstractMassiveMomentumGraph; kw...) =
  to_graphviz(to_graphviz_property_graph(g; kw...))

function to_graphviz_property_graph(g::AbstractMassiveMomentumGraph;
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
      e′ = add_edge!(pg, v, vertex(g, h), penwidth=string(mass(g, h) + 1))
    else
      e′ = add_edge!(pg, vertex(g, h), v, penwidth=string(mass(g, h) + 1))
    end
    set_eprops!(pg, e′, edge_label(g, edge_labels, h))
  end
  for (source, sink) ∈ zip(half_edge_pairs(g)...)
    (src, tgt) = (vertex(g, source), vertex(g, sink))
    e = add_edge!(pg, src, tgt, penwidth=string(mass(g, source) + 1))
    set_eprops!(pg, e, edge_label(g, edge_labels, e))
  end
  pg
end
















