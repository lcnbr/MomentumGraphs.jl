


@present SchMassiveMomentumGraph <: SchMomentumGraph begin
  Weight::AttrType
  mass::Attr(H, Vec)
end

@abstract_acset_type AbstractMassiveMomentumGraph <: AbstractMomentumGraph


@acset_type MassiveMomentumGraphGeneric(SchMomentumGraph, index=[:inv, :vertex, :sink,:momentum,:indep,:mass]) <: AbstractMassiveMomentumGraph

MassiveMomentumGraph=MassiveMomentumGraphGeneric{Real,Bool,FVector}

mass(g::AbstractMassiveMomentumGraph, args...) = subpart(g, args..., :mass)


function add_half_edge_pairs!(g::AbstractMassiveMomentumGraph, srcs::AbstractVector{Int},
  tgts::AbstractVector{Int}; indeps=falses(length(srcs)),mass=zeros(length(srcs)),kw...)

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
  set_subpart!(g, outs, :mass,mass)
  set_subpart!(g, ins, :mass,mass)

  first(outs):last(ins)
end

function add_dangling_edge!(g::AbstractMassiveMomentumGraph, v::Int; dir=true,momentum=FVector(value(variables(:p,length(dangling_edges(g))+1)[1])),indep=true,mass=0.0, kw...)

  H=add_part!(g, :H; vertex=v, inv=nparts(g,:H)+1,sink=dir,momentum=momentum,indep=indep,mass=mass,kw...)
end

function add_dangling_edges!(g::AbstractMassiveMomentumGraph, vs::AbstractVector{Int}; dirs::AbstractVector{Bool}=trues(length(vs)),momenta=FVector.(value.(variables(:p,(length(dangling_edges(g))+1):length(dangling_edges(g))+length(vs)))),indeps=trues(length(vs)),masses=zeros(length(vs)),kw...)
  neIn=length(dangling_edges(g))
  n, k = length(vs), nparts(g, :H)
  H=add_parts!(g, :H, n; vertex=vs, inv=(k+1):(k+n),sink=dirs, momentum=momenta,indep=indeps,mass=masses,kw...)
end

  

  









