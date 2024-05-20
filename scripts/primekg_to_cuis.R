devtools::load_all()
dirpath = system.file('extdata', package = 'kgraph')
df_kg = data.table::fread(file.path(dirpath, 'kg.csv'), data.table = FALSE)
df_kg_disease = subset(df_kg, relation == 'disease_disease')
# nrow(df_kg_disease)
# [1] 64388


fpath = file.path(dirpath, 'umls_mondo.csv')
df_mondo_cui = data.table::fread(fpath, data.table = FALSE)
# rm CN starting CUI codes, none in CUI to phecode mapping anyway
df_mondo_cui %<>% subset(!grepl('^CN', umls_id)) # only 1k out of 30k

df_merge = merge(df_kg_disease, df_mondo_cui, by.x = 'x_id', by.y = 'mondo_id')
df_merge %<>% merge(df_mondo_cui, by.x = 'y_id', by.y = 'mondo_id')
# nrow(df_merge)
# [1] 104874
df_cuis_pairs = df_merge
df_cuis_pairs %<>% subset(!duplicated(paste0(umls_id.x, umls_id.y)))
df_cuis_pairs %<>% subset(umls_id.x != umls_id.y)

df_cuis_pairs = df_cuis_pairs[c('umls_id.x', 'x_name', 'umls_id.y', 'y_name')]

dirpath = system.file('data', package = 'kgraph')
save(df_cuis_pairs, file = file.path(dirpath, 'df_cuis_pairs.rds'))

