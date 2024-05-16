fit_embeds_to_pairs = function(m_embeds, df_pairs,
  df_pairs_cols = c('PheCode.x', 'PheCode.y'),
  similarity = c('inprod', 'cosine', 'cov_simi', 'norm_inprod'),
  sparse_coding = c('none', 'epsilon'), add_concepts = NULL, ...) {

  similarity = match.arg(similarity)
  sparse_coding = match.arg(sparse_coding)

  # intersect pairs and embeds
  concepts_embeds = rownames(m_embeds)
  pairs_intersect = intersect(concepts_embeds,
                              unique(unlist(df_pairs[df_pairs_cols])))

  df_pairs %<>% subset(df_pairs[[df_pairs_cols[1]]] %in% pairs_intersect &
                       df_pairs[[df_pairs_cols[2]]] %in% pairs_intersect)

  # if not too much concepts, get similarity on all for projections
  if (nrow(m_embeds) > 4e3) {
      concepts_embeds = pairs_intersect
      if (!is.null(add_concepts)) concepts_embeds %<>% c(add_concepts)
  }

  m_simi = if (sparse_coding == 'epsilon') {
      stopifnot(dist_method %in% c('cosine', 'norm_inprod'))
      sparse_encode(m_embeds, dist_method = similarity, ...)
    } else {
      similarity_fun = switch(similarity, cosine = text2vec::cosine,
                              get(similarity))
      similarity_fun(m_embeds)
    }

  # we need as many false relations to test AUC,
  # make sure true pairs are also ordered
  df_pairs[df_pairs_cols] %<>% order_dataframe
  df_allpairs = df_pairs[df_pairs_cols] %>%
      rbind(gen_df_notpairs(pairs_intersect, .), .)

  # not very elegant right
  sims = apply(df_allpairs, 1, function(pairs) m_simi[pairs[[1]], pairs[[2]]])

  # fit roc
  truth = rep(0:1, each = nrow(df_pairs))
  roc_obj = pROC::roc(truth, sims, direction = '<', quiet = TRUE)

  # return what's needed for partial AUC and projections
  threshold_5fp = get_cutoff_threshold(roc_obj)
  df_projs = project_pairs(m_simi,
                           get_cutoff_threshold(roc_obj, 0.9))

  embeds_pairs_fit = list(roc = roc_obj, sims = sims, truth = truth,
                          threshold_5fp = threshold_5fp,
                          n_concepts = length(concepts_embeds),
                          df_projs = df_projs)

  embeds_pairs_fit
}

get_ppv = function(l_fit_embeds,
                   threshold = get_cutoff_threshold(l_fit_embeds$roc)) {

  preds = ifelse(l_fit_embeds$sims > threshold, 1, 0)
  sum(preds == 1 & preds == l_fit_embeds$truth) / sum(preds == 1)
}

get_cutoff_threshold = function(roc_obj, specificity_lvl = 0.95) {

  idx = abs(roc_obj$specificities - specificity_lvl) %>% which.min
  roc_obj$thresholds[idx]
}

cov_simi = function(m_data) cov(t(m_data))

get_known_pairs = function() {

  dirpath = system.file('evaluations', package = 'psychclust')
  df_pairs = get(load(file.path(dirpath, 'pairs_arranged.Rdata')))
  df_pairs = subset(df_pairs, group == 'cui-cui')[c('code1', 'code2')]
}

project_pairs = function(m_simi, threshold) {

  simi_dimnames = dimnames(m_simi)
  m_simi = ifelse(lower.tri(m_simi), m_simi, 0)
  dimnames(m_simi) = simi_dimnames
  df_preds = reshape2::melt(m_simi) %>% subset(value > threshold)

  df_preds = df_preds[order(df_preds$value, decreasing = TRUE), ]

  setNames(df_preds, c('concept1', 'concept2', 'weight'))
}

setdiff_dataframe = function(df_x, df_y, cols = 1:2) {
  df_x %>% subset(!Reduce('paste0', .[cols]) %in% Reduce('paste0', df_y[cols]))
}

order_dataframe = function(df_x, cols = 1:2, relevant_pattern = NULL) {

  # put strings matching relevant patterns in first column
  if (!is.null(relevant_pattern)) {

    df_x %<>% subset(grepl(relevant_pattern, .[[cols[1]]]) |
                     grepl(relevant_pattern, .[[cols[2]]]))

    uniq_lvls = unique(unlist(df_x[cols]))
    relevant_lvls = grep(relevant_pattern, uniq_lvls)
    uniq_lvls = c(uniq_lvls[relevant_lvls], uniq_lvls[-relevant_lvls])

    df_x[cols] %<>% lapply(function(lvls) factor(lvls, uniq_lvls) %>%
                           as.numeric)
    df_x %<>% order_dataframe
    df_x[cols] %<>% lapply(function(lvls) uniq_lvls[lvls])

  } else {
    df_x[cols] = apply(df_x[cols], 1, sort) %>%
      t %>% as.data.frame %>% setNames(names(df_x[cols]))
  }

  df_x
}

gen_df_notpairs = function(ids, df_pairs) {

  set.seed(1)
  df_notpairs = sample(ids, nrow(df_pairs) * 4, replace = TRUE) %>%
      matrix(ncol = 2) %>% as.data.frame %>% setNames(names(df_pairs[1:2]))

  df_notpairs %>% subset(.[[1]] != .[[2]]) %>% order_dataframe %>%
    setdiff_dataframe(df_pairs) %>% head(nrow(df_pairs))
}
