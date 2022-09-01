module TestKarnak
  using DirectedHalfEdgeGraphs
  using MomentumGraphs
  using Karnak
  using Colors
  using Catlab.CategoricalAlgebra.CSets
  using Catlab.Graphs
  import Graphs: AbstractGraph
  Graphs.AbstractGraph=BasicGraphs.AbstractGraph
  g = Graph()
add_vertices!(g, 4)
add_edge!(g, 1, 2)
add_edge!(g, 1, 3)
add_edge!(g, 2, 3)
add_edge!(g, 1, 4)



@drawsvg begin
    background("grey10")
    sethue("orange")
    drawgraph(g, vertexlabels = [1, 2, 3, 4])
end 600 300
end