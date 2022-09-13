module Fields
using ...Form
export AbstractField,AbstractScalarField,ScalarField,AbstractVectorField,VectorField,AbstractTensorField,TensorField,AbstractGaugeVectorField,GaugeVectorField,decoration,ismassless,AdjointField,basefield,color,isadjoint,isselfadjoint

abstract type AbstractField end

(::Type{T})(f::T) where {T<: AbstractField} =f

struct AdjointField{S} <: AbstractField
  parent::S
  AdjointField{S}(parent::S) where S<:AbstractField = new{S}(parent)
end

(::Type{S})(f::AdjointField{S}) where {S<:AbstractField}=f.parent
AdjointField(parent::S) where S<:AbstractField = AdjointField{S}(parent)
AdjointField{S}(arg...) where S<: AbstractField = AdjointField(S(arg...))

Base.show(io::IO, ::MIME"text/plain", ::Type{AdjointField{Kind}}) where Kind = print(io,(string("Adjoint ",Kind)))

Base.show(io::IO, ::MIME"text/plain", f::AdjointField{Kind}) where Kind = print(io,(string("Adjoint ",f.parent)))

Base.show(io::IO, f::AdjointField{Kind}) where Kind = print(io,(string("Adjoint ",f.parent)))

decoration(::AbstractField) = :none

Base.adjoint(f::T) where T<: AbstractField = Base.adjoint(T)(f)
Base.adjoint(::Type{AdjointField{S}}) where S  = S
Base.adjoint(::Type{T}) where T<: AbstractField = AdjointField{T}

isselfadjoint(::Type{T}) where T<: AbstractField = T==adjoint(T)
isselfadjoint(::T) where T<: AbstractField = isselfadjoint(T)
isadjoint(::T) where T<: AbstractField = isadjoint(T)
isadjoint(::Type{T}) where T<: AbstractField = isselfadjoint(T)
isadjoint(::Type{T}) where T<: AdjointField = true

basefield(f::T) where T<: AbstractField = f
basefield(f::T) where T<: AdjointField = f.parent
basefield(::Type{T}) where T<: AbstractField = T
basefield(::Type{AdjointField{S}}) where S = S

color(::T) where T<: AbstractField = color(T)
color(::Type{AdjointField{S}}) where S = color(S)
color(::Type{T}) where T<: AbstractField = "black"

mass(::T) where T<: AbstractField = mass(T)
mass(::Type{AdjointField{S}}) where S = mass(S)
mass(::Type{T}) where T<: AbstractField = :none

ismassless(::T) where T<: AbstractField = ismassless(T)
ismassless(::Type{AdjointField{S}}) where S = ismassless(S)
function ismassless(::Type{T}) where T<: AbstractField 
  @show mass(T)
   ((@show mass(@show T)) === :none)
end
abstract type AbstractScalarField <: AbstractField end
struct ScalarField{Kind} <: AbstractScalarField 
  name::Symbol
  ScalarField{Kind}(name::Symbol=Kind) where Kind = new(name)
end



Kind = ScalarField{Kind} where Kind

ScalarField(name::Symbol) = ScalarField{name}()
(::Type{T})(ID::Int,name...) where T<: AbstractScalarField= T(name...)

(::Type{T})(f::AbstractScalarField) where T <: AbstractScalarField=T(f.name)

Base.show(io::IO, ::MIME"text/plain", ::Type{ScalarField{Kind}}) where Kind = print(io,(string("Scalar: ",Kind)))

Base.show(io::IO, ::Type{ScalarField{Kind}}) where Kind = print(io,(string(Kind)))

function Base.show(io::IO,f::ScalarField{Kind}) where Kind
  if Kind===f.name
    print(io,(string(Kind)))
  else
    print(io,(string(Kind," ",f.name)))
  end
end



abstract type AbstractVectorField <: AbstractField end
struct VectorField{Kind} <: AbstractVectorField 
  name::Symbol      
  index::FIndex
  VectorField{Kind}(index::FIndex,name::Symbol=Kind) where Kind =new(name,index)
end

VectorField(index::FIndex,name::Symbol) = VectorField{name}(index) 
(::Type{T})(ID::Int,name...) where T<: AbstractVectorField= T(FIndex("mu",ID),name...)

(::Type{T})(f::VectorField) where T <: AbstractVectorField=T(f.index,f.name)

Base.show(io::IO, ::MIME"text/plain", ::Type{VectorField{Kind}}) where Kind = print(io,(string("Vector: ",Kind)))

Base.show(io::IO, ::Type{VectorField{Kind}}) where Kind = print(io,(string(Kind)))

function Base.show(io::IO,f::VectorField{Kind}) where Kind
  if Kind===f.name
    print(io,(string(Kind)),"(",f.index,")")
  else
    print(io,(string(Kind," ",f.name,"(",f.index,")")))
  end
end

abstract type AbstractTensorField <: AbstractField end
struct TensorField{Kind} <: AbstractTensorField 
  name::Symbol
  indices::AbstractVector{FIndex}
  TensorField{Kind}(indices::AbstractVector{FIndex},name::Symbol=Kind) where Kind =new(name,indices)
end 

TensorField(indices::AbstractVector{FIndex},name::Symbol) = TensorField{name}(indices) 
(::Type{T})(f::AbstractTensorField) where T <: AbstractTensorField =T(f.indices,f.name)
(::Type{T})(ID::Int,size::Int,name...) where T <: AbstractTensorField =T( FIndex.([Symbol(string(greek(i),subscript(ID))) for i in 1:size]),name...)

greek(i) = Char(0x03b1 + (i-1)%25)*(i>25 ? superscript((i÷25)-1) : "")
subscript(i) = join(Char(0x2080 + d) for d in reverse!(digits(i)))
superscript(i) = join(Char(0x2070 + d) for d in reverse!(digits(i)))

Base.show(io::IO, ::MIME"text/plain", ::Type{TensorField{Kind}}) where Kind = print(io,(string("Tensor: ",Kind)))

abstract type AbstractGaugeVectorField <: AbstractVectorField end
struct GaugeVectorField{Kind} <: AbstractGaugeVectorField
  name::Symbol
  index::FIndex
  function GaugeVectorField{Kind}(index::FIndex,name::Symbol=Kind) where Kind
    @info "construct"
    new(name,index)
  end
end  
(::Type{T})(f::AbstractGaugeVectorField) where T <: AbstractGaugeVectorField=T(f.index,f.name)
GaugeVectorField(index::FIndex,name::Symbol) = GaugeVectorField{name}(index)
function Base.show(io::IO,f::GaugeVectorField{Kind}) where Kind
  if Kind===f.name
    print(io,(string(Kind)),"(",f.index,")")
  else
    print(io,(string(Kind," ",f.name,"(",f.index,")")))
  end
end

Base.show(io::IO, ::MIME"text/plain", ::Type{GaugeVectorField{Kind}}) where Kind = print(io,(string("Gauge: ",Kind)))
  
end