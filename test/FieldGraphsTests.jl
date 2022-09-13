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
using MomentumGraphs.Fields
import MomentumGraphs.Fields: mass,color

FieldGraph()

Photon=GaugeVectorField{:Photon}
Base.adjoint(::Type{Photon})=Photon
Phi1 =  ScalarField{:Phi1}
Phi1c = adjoint(Phi1)
Phi2 = ScalarField{:Phi2}
Phi2c=  adjoint(Phi2)

Fields.mass(::Type{Phi1})=:m1
Fields.mass(::Type{Phi2})=:m2
color(::Type{Photon})="grey"
color(::Type{Phi1})="blue"
color(::Type{Phi2})="red"
color(::Type{Phi1c})="deepskyblue"
color(::Type{Phi2c})="lightcoral"

fg = FieldGraph([-2,1,-1,-4,3,-3,8,2,5,6,4,7],
[0,2,0,0,4,0,7,1,6,5,3,8],
[1,1,1,2,2,2,3,3,3,4,4,4],
[Phi1c,Photon,Phi1,Phi2c,Photon,Phi2,Phi2c,Photon,Phi2,Phi2c,Photon,Phi2];strict=true)
to_graphviz(fg)
map.(momentum(fg),map(x->x.index,field(fg)))
fg = FieldGraph(5)

@test add_vertex!(fg) == 6
#add_edges!(fg, [1,2,4,3], [2,4,3,1])
@test add_edges!(fg, [1, 1, 1, 2, 2, 4, 3, 5, 6], [2, 5, 3, 5, 4, 6, 4, 6, 5]) == 1:18
@test add_dangling_edges!(fg, [1, 2, 3, 4], dirs=[true, true, false, false]) == 19:22

fg
to_graphviz(fg)


fg = FieldGraph(5)
@test add_vertex!(fg) == 6

γ=VectorField{:γ}
Base.adjoint(::Type{γ})=γ
ψ =  VectorField{:ψ}
ψ̄ = adjoint(ψ)
ϕ=GaugeVectorField{:ϕ}
ϕ̄ = adjoint(ϕ)

Fields.mass(::Type{ϕ})=:m1
Fields.mass(::Type{ψ})=:m2

color(::Type{γ})="grey"
color(::Type{ψ})="blue"
color(::Type{ϕ})="red"



fieldtypes=[ϕ,ψ,γ,ϕ,ψ,γ,ϕ,ψ,γ]


#add_edges!(fg, [1,2,4,3], [2,4,3,1])

@test add_edges!(fg, [1, 1, 1, 2, 2, 4, 3, 5, 6], [2, 5, 3, 5, 4, 6, 4, 6, 5];fieldtypes) == 1:18
@test add_dangling_edges!(fg, [1, 2, 3, 4]; dirs=[true, true, false, false],fieldtypes=[ϕ,ϕ,ψ,ψ]) == 19:22

fg
to_graphviz(fg)
diag1=(
  H           = [-2,1,3,-1,-4,2,4,-3],
  Hdual       = [0,2,4,0,0,1,3,0],
  vertex      = [1,1,1,1,2,2,2,2],
  field       = [:phi1c,:photon,:photon,:phi1,:phi2c,:photon,:photon,:phi2]
  )

FieldGraph([-2,1,3,-1,-4,2,4,-3],[0,2,4,0,0,1,3,0],[1,1,1,1,2,2,2,2],[:phi1c,:photon,:photon,:phi1,:phi2c,:photon,:photon,:phi2];dualDict,massDict)

spanningTree = subtree(fg, dfs_parents(fg, 4, all_neighbors, by=(x -> Base.sort(x, by=e -> mass(fg, e), rev=true))))
to_graphviz(spanningTree, node_labels=true, edge_labels=true; field_colors)

fg = Subobject(fg, H=half_edges(fg), V=vertices(fg))
to_graphviz(sfg \ spanningTree;field_colors)
sfg \ spanningTree |> hom |> dom


momentum_equations(fg)
momentum_equations_solved(fg;by=(x -> Base.sort(x, by=e -> mass(fg, e), rev=false)))


using Catlab.Graphs.GraphAlgorithms
connected_component_projection_bfs(fg)
connected_component_projection(fg)

using CSetAutomorphisms

nickel_index(call_nauty(fg).cset)

import MomentumGraphs.FieldGraphs: GaugeVectorField, decoration,dual,ScalarField,mass,ismassless,AbstractField

abstract type AbstractField end
abstract type AbstractScalarField <: AbstractField end


using MomentumGraphs.FieldGraphs: AbstractVectorField,VectorField,ScalarField,mass,ismassless,AbstractField,TensorField
TensorField{:test}(5,26)
TensorField(5,26,:test)

VectorField{:s}(3)
ScalarField(:s)

function VectorField{Kind}(name,ID::Int) where Kind
  ind=FIndex("mu",ID)
  VectorField{Kind}(name,ind)
end
VectorField{:Psi2}(VectorField{:Psi}(:s,2))
FIndex("mu",2)
struct ScalarField{Kind} <: AbstractScalarField 
  name::Symbol
end

#(::Type{T})(f::F) where {T<: ScalarField,F <: ScalarField} =T(f.name) 
Base.convert(::Type{T}, a::F) where {T<: ScalarField,F <: ScalarField} = T(a.name) 
decoration(::AbstractField) = :none
(::Type{ScalarField{Kind}})(f::ScalarField) where Kind =ScalarField{Kind}(f.name)

Base.show(io::IO, ::MIME"text/plain", ::Type{ScalarField{Kind}}) where Kind = print(io,(string("Scalar: ",Kind)))




Base.adjoint(f::T) where T<: AbstractField = Base.adjoint(T)(f)
Base.adjoint(::Type{T}) where T<: AbstractField = T


Phi1 = ScalarField{:Phi1}
Phi1c = ScalarField{:Phi1c}
Phi2 = ScalarField{:Phi2}
Phi2c=  ScalarField{:Phi2c}
Phi3 = ScalarField{:Phi3}
Phi1("s")
a=Phi1(Phi1(:ψ))
(a)'



decoration(::Photon)=:snake
Base.adjoint(::Type{Phi1}) = Phi1c
dual(::Type{Phi1c}) = Phi1
dual(::Type{Phi2}) = Phi2c
dual(::Type{Phi2c}) = Phi2
dual(::Type{Photon}) = Photon


mass(::Type{Phi1})=:m1
mass(::Type{Phi1c})=:m1
mass(::Type{Phi2})=:m2
mass(::Type{Phi2c})=:m2
mass(::Type{Photon})=:none


ismassless(Photon)
ismassless(Phi1)
ismassless(Phi1c) 
ismassless(::Type{Phi2})=false
ismassless(::Type{Phi2c})=false


end