build_kgraph = function(selected_concept, df_weights, df_dict,
                        rm_single_groups = TRUE, add_drugs = FALSE,
                        display_val_digits = 3, str_other = 'Other') {

  message(paste('Building', selected_concept))

  df_weights %<>% subset(concept1 == selected_concept |
                         concept2 == selected_concept)
  df_weights %<>% order_dataframe(relevant_pattern = selected_concept)

  df_merge = merge(df_weights, df_dict, by.x = 'concept2', by.y = 'id')

  # link all nodes to central node
  first_order_links = df_merge[c('concept1', 'concept2', 'weight')] %>%
      setNames(c('from', 'to', 'weight'))

  # link nodes to groups
  if (!is.null(df_merge$group)) {
    l_df_merge = merge_per_groups(df_merge, rm_single_groups, str_other)
    second_order_links = switch_clusters_type(l_df_merge, selected_concept,
                                              type = 2) 
  } else {
    l_df_merge = list(df_merge)
    second_order_links = NULL
  }

  # link supplementary nodes
  if (add_drugs) {

    df_drugs = subset(df_merge, !is.na(drug_name))
    l_df_drugs = split(df_drugs, df_drugs$drug_name)

    third_order_links = switch_clusters_type(l_df_drugs, selected_concept,
                                             type = 3, from_col = 'drug_name') 

    third_order_groups = switch_clusters_type(l_df_drugs, selected_concept,
                                              type = 3, from_col = 'drug_name',
                                              to_col = 'group')

    # in case of rm_single_groups
    fix_rm_single_grps = !third_order_groups$to %in% second_order_links$from
    third_order_groups$to[fix_rm_single_grps] = str_other
    third_order_links %<>% rbind(third_order_groups) %>% unique

    wts = diff(range(df_merge$weight)) / c(5, 4, 3, 2)
    phase_idx = match(third_order_links$from, df_drugs$drug_name)
    drug_phase = df_drugs$max_phase[phase_idx]
    third_order_links$weight = wts[drug_phase]

  } else {
    third_order_links = NULL
  }

  df_nodes = build_kgraph_nodes(first_order_links, second_order_links,
                                third_order_links, df_weights, df_dict,
                                l_df_merge, selected_concept,
                                display_val_digits) 

  df_links = rbind(first_order_links, second_order_links, third_order_links)

  stopifnot(nrow(df_nodes) == length(unique(unlist(df_links[1:2]))))

  list(df_links = df_links, df_nodes = df_nodes)
}

build_kgraph_nodes = function(first_order_links, second_order_links,
                              third_order_links, df_weights, df_dict,
                              l_df_merge, selected_concept,
                              display_val_digits) {

  df_root = df_dict[c('id', 'desc', 'color')] %>%
      subset(id %in% selected_concept)

  df_merge = do.call(rbind, l_df_merge)
  leaf_idxs = match(first_order_links$to, df_merge$concept2)
  df_leafs = df_merge[leaf_idxs, c('concept2', 'desc', 'color')]
  colnames(df_leafs)[1] <- "id"
  
  if (!is.null(second_order_links)) {
    group_idxs = match(unique(second_order_links$from), df_dict$group)
    df_groups = unique(second_order_links$from) %>%
        cbind(id = ., desc = ., color = df_dict$color[group_idxs])
    df_groups$color = 'Groups'
  } else {
    df_groups = NULL
  }
  
  if (!is.null(third_order_links)) {

    uniq_drugs = unique(third_order_links$from)

    drugs_idxs = match(uniq_drugs, df_merge$drug_name)
    max_phase = df_merge$max_phase[drugs_idxs]

    df_drugs = cbind(id = uniq_drugs, desc = max_phase,
              color = df_merge$color[drugs_idxs])

    df_drugs$color = 'Drugs'

  } else {
    df_drugs = NULL
  }

  df_nodes = rbind(df_root, df_leafs, df_groups, df_drugs) %>%
      setNames(c('word', 'desc', 'clusters'))

  display_val_idx = match(df_nodes$word, df_weights$concept2)
  df_nodes$display_val <- round(df_weights$weight[display_val_idx],
                                display_val_digits)
  df_nodes$display_val[1] = NA

  wts = first_order_links$weight
  
  if (!is.null(df_groups)) {
    wts %<>% c(second_order_links$weight[!duplicated(second_order_links$from)])
  }
  
  if (!is.null(df_drugs)) {
    wts %<>% c(third_order_links$weight[!duplicated(third_order_links$from)])
  }

  # add selected_concept as first weight
  df_nodes$weight = wts %>% c(max(.) + (diff(range(.)) / 4), .)
  df_nodes$selected_concept = selected_concept

  ## reduce group edges
  #group_edges = df_links$from %in% subset(df_nodes, clusters == 'Groups')$word
  #df_links$weight[group_edges] %<>% `/`(2)

  df_nodes
}

# three different ways of organizing groups
# NOTE could improve type 3 by removing group-target links when multiple
## targets matched by children
# but type 2 probably better by default
# could improve type 2 by shortening outliers
switch_clusters_type = function(l_df_merge, selected_concept, type = 2,
                                from_col = 'group', to_col = 'concept2') {

  # either get same group weighting for all group nodes
  if (type == 1) {

    df_merge = do.call(rbind, l_df_merge)
    second_order_links = df_merge[c(from_col, to_col)] %>%
        setNames(c('from', 'to'))

    second_order_links$weight = median(df_merge$weight)

  # or get group node weights per group
  } else {

    second_order_links = lapply(l_df_merge, links_by_cluster, to_col) %>%
        reshape2::melt(id = to_col) %>%
        setNames(c('to', 'var', 'weight', 'from'))

    second_order_links = second_order_links[c(4, 1, 3)]

    # add an edge from group to central node (useful for small distant groups)
    if (type == 3) {

      third_order_links = second_order_links
      third_order_links$to = selected_concept
      third_order_links %<>% unique
      second_order_links %<>% rbind(third_order_links)
    }
  }

  second_order_links
}

links_by_cluster = function(df_merge_grp, to_col = 'concept2') {
   data.frame(df_merge_grp[[to_col]], max(df_merge_grp$weight)) %>%
   # data.frame(df_merge_grp[[to_col]],
   #    diff(range(df_merge_grp$weight)) / 1.5) %>%
      setNames(c(to_col, 'weight'))
}

# NOTE to sigmagraph
get_legend = function(colors_mapping, clusters) {

  colors_mapping %<>% subset(group %in% clusters)
  colors_mapping %<>% cbind(data.frame(x = 1, y = 1))
  colors_mapping$group %<>% factor(unique(.))

  gglegend = ggplot2::ggplot(colors_mapping,
                             ggplot2::aes(x, y, color = group)) +
      ggplot2::geom_point(size = 10) +
      ggplot2::scale_color_manual(name = NULL, values = colors_mapping$color) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.text.position = 'top',
                     legend.title = ggplot2::element_text(size = 20),
                     legend.text = ggplot2::element_text(size = 15))

  gglegend = cowplot::get_legend(gglegend)
}

get_color_map = function(color_levels) {

  palette = RColorBrewer::brewer.pal(length(color_levels), 'Paired')

  sigmagraph:::get_color_map(color_levels, palette)
}
