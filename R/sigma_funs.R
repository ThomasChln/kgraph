multiline_labels = function(df_nodes, display_val_str = '\nP-value: ',
                            replace_codes = TRUE) {

  if (!replace_codes) {
    df_nodes$label = df_nodes$desc %>% {
        ifelse(df_nodes$word != .,
               paste0(ifelse(grepl('^[0-9]$', .), 'Max. phase: ', 'Label: '),
                      ., '\n'),
               '')
    }

    df_nodes$label %<>% paste0('Group: ', df_nodes$clusters)

  } else {

    df_nodes$label = df_nodes$desc %>% ifelse(is.na(.), df_nodes$word, .)
  }

  val_labels = df_nodes$display_val %>%
      { ifelse(!is.na(.), paste0(display_val_str, df_nodes$display_val), '') }

  df_nodes$label %<>% paste0(val_labels)

  df_nodes
}

get_kgraph = function(selected_concepts, df_weights, df_dict) {

  kgraphs = selected_concepts %>% setNames(., .) %>%
    lapply(build_kgraph, df_weights, df_dict)

  kgraph = merge_kgraphs(kgraphs, spring_weights = TRUE, df_dict,
                         display_val_str = '\nCosine similarity: ')

  kgraph$df_links = kgraph %$%
	  sigmagraph:::highlight_multiple_connected(df_links, selected_concepts)

  kgraph$df_nodes %<>% multiline_labels('\nCosine similarity: ')

  kgraph
}

get_sgraph = function(l_graph, colors_mapping, label_attrs = 'label') {

  igraph = sigmagraph:::l_graph_to_igraph(l_graph)

  sgraph = sigmagraph:::sigma_mutuals(igraph, niter = 1e3,
				                      node_size = 'weight',
                                      label = label_attrs,
				                      color_map = colors_mapping,
                                      layout = igraph::layout_with_kk(igraph))
  
  sgraph %<>% sigmagraph::add_edge_color(one_color = l_graph$df_links$color)
  sgraph %<>% sigmagraph:::add_edge_zindex(zindex = l_graph$df_links$zindex)
}
