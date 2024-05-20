m_embeds = get(load('data/fusion_glove_fit_nile.rds'))

#df_dict = get(load('data/df_cui_mapping_nile.rds'))

#df_dict = df_dict[na.omit(match(rownames(m_embeds), df_dict$id)), ]
#df_dict = df_dict[c('id', 'label')]
#save(df_dict, file = 'data/df_dict.rds')

fpath = 'inst/portable_NILE/project/test/input/test_dict.txt'
df_dict = data.table::fread(fpath, data.table = FALSE)
df_dict = df_dict[c(2, 1, 4, 5)] %>%
       	setNames(c('id', 'label', 'group', 'group_details'))

map_idxs = match(rownames(m_embeds), df_dict$id)
df_dict = df_dict[na.omit(map_idxs), ]
m_embeds = m_embeds[!is.na(map_idxs), ]

save(df_dict, file = 'data/df_dict.rds')
save(m_embeds, file = 'data/m_embeds.rds')

