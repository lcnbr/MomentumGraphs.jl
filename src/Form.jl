module Form
using SymbolicUtils
using Symbolics
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