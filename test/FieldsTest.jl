module TestFields

using MomentumGraphs.Form

using MomentumGraphs.Fields
TensorField{:test}(5,26)
TensorField(5,26,:test)

VectorField{:s}(3)
ScalarField(:s)

const Phi1    = ScalarField{:Phi1}
Photon=GaugeVectorField{:Photon}
Base.adjoint(::Type{Photon})=Photon
p =Photon()
adjoint(p)
Phi1c = AdjointField{ScalarField{:Phi1}}
adjoint(Phi1c)
AdjointField{ScalarField{:Phi1}}(1)
adjoint(Phi1(1,:test))
Photon(1)
endswith
field(Phi1c(1))
field()
Phi2 = ScalarField{:Phi2}
Phi2c=  ScalarField{:Phi2c}
Phi3 = ScalarField{:Phi3}
Phi1(:S)
a=Phi1(Phi1(:ψ))
adjoint.([a,a,a])
(a)'


struct test
  fielt
end
test(test(4))

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


abstract type A end
abstract type B <: A end  
struct S<:B
  field1
  field2
end
struct R<:B
  field1
  field2
end

function (::Type{T})(f::T) where {T<: A} 
  @info "no copy"
  f
end

function (::Type{T})(f::B) where T <: B
  @info  "copy"
  T(f.field1,f.field2)  
end

r=R(2,3)
S(r)
R(r)

end