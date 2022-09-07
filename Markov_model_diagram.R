library(DiagrammeR)
library(diagram)

a_graph <-
  create_graph() %>%
  add_node(label = "Asymptomatic",
           node_aes = node_aes(fontsize = 5)) %>% # 1
  add_node(label = "Progressive \ndisease",
           from = 1,
           node_aes = node_aes(fontsize = 6),
           edge_aes = edge_aes(label = "tpProg*(1 - effect)", fontsize = 6, len = 5)) %>%
  add_node(label = "Dead",
           from = c(1,2),
           node_aes = node_aes(fontsize = 8),
           edge_aes = edge_aes(label = c("tpDn", "tpDcm + tpDn"), fontsize = 6, len = 5)) %>%
  add_edge(from = 1, to = 1, edge_aes = edge_aes(label = " ")) %>%
  add_edge(from = 2, to = 2,
           edge_aes = edge_aes(fontsize = 1)) %>%
  add_edge(from = 3, to = 3,
           edge_aes = edge_aes(fontsize = 10))

render_graph(a_graph, layout = "nicely")
render_graph(a_graph, layout = "circle")


# library(DiagrammeRsvg)
export_graph(graph = a_graph,
             file_name = "markov_model_diagram.png",
             file_type = "PNG")




