module test

using Latexify
import Latexify: unicodedict,unicode2latex

using Test

  using MomentumGraphs.Form

  FIndex()
  FIndex(:μ₁)
  FIndex("μ",1).ID
  FIndex("mu",1).ID

  greek(i) = Char(0x03b0 + i%26)
  subscript(i) = join(Char(0x2080 + d) for d in reverse!(digits(i)))

  FIndex.([Symbol(string(greek(i),subscript(333))) for i in 1:25])
  greek(25)
end