
options(java.parameters=c("-XX:+UseConcMarkSweepGC", "-Xmx4096m"))


library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(diffr)
library(purrr)
library(xlsx)
library(seqinr)
library(ggridges)


random_seed <- 1
low_percentile <- 0.002
high_percentile <- 0.025
volcano_threshold_linetype <- 2
volcano_threshold_linesize <- 0.5
enrich_c <- .05 # increase to make it more "round" or decrease to make it more "edgy" :-)
enrich_s0 <- log2(2) # two fold change
enrich_pvalue <- 0.05

exp_name_level <- c('E. coli', 'B. subtilis', 'H. salinarum', 'T. brucei', 'T. thermophila',
                    'S. cerevisiae', 'S. pombe', 'C. elegans', 'C. elegans 1', 'C. elegans 2', 
                    'H. sapiens (HeLa)', 
                    'H. sapiens (HEK293)',  'H. sapiens (HeLa) ENSEMBL', 
                    'H. sapiens (HEK293) ENSEMBL', 'Z. mays', 'B. oleracea', 'A. thaliana')

out_folder <- 'specific_databases'
if(!dir.exists(out_folder)) dir.create(out_folder)

datasets_df <- read.delim('datasets_df_spec.txt', stringsAsFactors=FALSE) %>%
  as_tibble()

all_pg <- 
  read.delim('all_pg_spec_incl_ensembl_majority.txt', stringsAsFactors=FALSE) %>% 
  filter(exp.name != 'C. elegans 1') %>% 
  left_join(datasets_df %>% select(exp.name=species, short.name)) %>% 
  mutate(exp.name=if_else(exp.name == 'C. elegans 2', 'C. elegans', exp.name)) %>% 
  as_tibble()

all_pg %>% 
  select(Protein.IDs, short.name) %>% 
  mutate(Protein.IDs=strsplit(Protein.IDs, ';')) %>% 
  unnest(Protein.IDs) %>% 
  distinct_all() %>% 
  group_split(short.name) ->
  detected_ids_split
if(!dir.exists('individual_species_ids')) {
  dir.create('individual_species_ids')
  for(df in detected_ids_split) {
    writeLines(
      df$Protein.IDs, 
      file.path('individual_species_ids', paste0(unique(df$short.name),'.txt'))
    )    
  }
}

orthoMCL <- 
  read.delim('groups_OrthoMCL-CURRENT_table.txt', stringsAsFactors=FALSE) %>% 
  mutate(ortho_Protein.IDs=Protein.IDs,
         Protein.IDs=sub('^....\\|','', Protein.IDs)) %>% 
  as_tibble()
gene_names <- read.delim('individual_species_ids/gene_names.txt', stringsAsFactors=FALSE) %>% 
  as_tibble()
gene_name_order <- c("hsap", "scer", "spom", "cele", "ecol", "atha", "bsub", "halo", 
                     "tbrt", "tetr", "zmay")


orthoMCL %>% 
  mutate(species=sub('([^\\|]+)\\|.*', '\\1', ortho_Protein.IDs)) %>% 
  left_join(full_gene_info) %>%
  # distinct() %>% 
  mutate(species=factor(species, levels=gene_name_order)) %>% 
  arrange(species) %>% 
  group_by(ortho_group) %>% 
  slice(1) %>% 
  mutate(Gene.names_long=if_else(is.na(Gene.names), 
                                 as.character(NA), 
                                 paste0(species, '_', Gene.names))) %>% 
  filter(!is.na(Gene.names_long)) %>% 
  select(ortho_group, Gene.names_long) ->
  ortho_gene_names

orthoMCL %>% 
  left_join(all_pg %>% 
              mutate(id=paste(exp.name, id)) %>%
              select(Protein.IDs, id) %>%
              mutate(Protein.IDs=strsplit(Protein.IDs, ';')) %>%
              unnest(Protein.IDs)) %>% 
  group_by(ortho_group) %>% 
  summarise(in_ortho_group=n(),
            detected_proteins=sum(!is.na(id)),
            detected_proteinGroups=sum(!is.na(unique(id))),
            detected_pg_ids=paste0(unique(id[!is.na(id)]), collapse=';'),
            all_ortho_ids=paste0(ortho_Protein.IDs, collapse=';'),
            detected_ortho_ids=paste0(ortho_Protein.IDs[!is.na(id)], collapse=';')) %>% 
  write.table(file.path(out_folder, 'detected_ortholog_groups.txt'),
              row.names=FALSE, sep='\t')

ortho_summary <- 
  read.delim('AllGroups_Summary.txt', stringsAsFactors=FALSE) %>% 
  as_tibble() %>%
  mutate(Keywords=if_else(Keywords == 'N/A', as.character(NA), Keywords)) %>%
  select(-X)

all_pg %>%
  select(exp.name, Protein.IDs, id) %>%
  mutate(Protein.IDs=strsplit(Protein.IDs, ';')) %>%
  unnest(Protein.IDs) %>%
  left_join(orthoMCL, by='Protein.IDs') %>%
  left_join(ortho_summary, by=c('ortho_group'='Ortholog.Group')) %>% 
  group_by(exp.name, id) %>% 
  summarise(orthoMCL=paste0(na.omit(unique(ortho_group)), collapse=';'),
            ortho_keywords=paste0(na.omit(unique(Keywords)), collapse='; '),
            ortho_Protein.IDs=paste0(na.omit(unique(ortho_Protein.IDs)), collapse=';')) ->
  ortho_db
  

all_pg %>%
  filter(Reverse != '+',
         Potential.contaminant != '+',
         Only.identified.by.site != '+',
         Peptides >= 2,
         Unique.peptides >= 1,
         Razor...unique.peptides >= 2) %>%
  left_join(ortho_db, by=c('exp.name','id')) %>% 
  mutate(exp.name=factor(exp.name, levels=exp_name_level))->
  pg

cat(sprintf('There are %d Proteins with missing ortholog group!\n',
            sum(!grepl('^OG6.*', pg$orthoMCL))))

jnk <- read.delim('pg_all_orthoMCL.txt')
jnk %>% 
  filter(pulldown.date != 20201103 | exp.name != 'C. elegans') %>% 
  select(exp.name,id) %>% 
  mutate(dataset='orthoMCL') %>% 
  bind_rows(pg %>% 
              select(exp.name,id) %>% 
              mutate(dataset='specific')) %>% 
  group_by(dataset, exp.name) %>% 
  summarise(count=n()) %>% 
  ggplot(aes(exp.name, count, fill=dataset)) +
  geom_col(position=position_dodge(width=.9)) + 
  geom_text(aes(label=count), vjust=0, position=position_dodge(.9),
            size=3) +
  labs(title='protein group count after filtering (2peps, 1unique)',
       subtitle='comparing specific or orthoMCL provided databases') +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(file.path(out_folder, 'protein_count_specific_vs_orthoMCL.pdf'),
       width=11.69, height=8.27)

jnk %>% 
  filter(pulldown.date != 20201103 | exp.name != 'C. elegans') %>% 
  select(exp.name,id) %>% 
  mutate(dataset='orthoMCL') %>% 
  bind_rows(pg %>% 
              select(exp.name,id) %>% 
              mutate(dataset='specific')) %>% 
  group_by(dataset, exp.name) %>% 
  summarise(count=n()) %>% 
  pivot_wider(names_from=dataset, values_from=count) %>% 
  mutate(`specific - orthoMCL`=specific - orthoMCL) %>% 
  ggplot(aes(exp.name, `specific - orthoMCL`)) +
  geom_col() + 
  geom_text_repel(aes(label=`specific - orthoMCL`), size=3,direction='y') +
  labs(title='protein group count difference (2peps, 1unique)',
       subtitle='subtracting orthoMCL from specific database counts') +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(file.path(out_folder, 'protein_count_difference_specific_vs_orthoMCL.pdf'),
       width=11.69, height=8.27)

pg %>% 
  group_by(exp.name) %>% 
  summarise(n=sum(!grepl('^OG6.*', orthoMCL))) %>% 
  ggplot(aes(exp.name, n)) +
  geom_col() +
  labs(title='Protein group count for missing ortholog groups') +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave(file.path(out_folder, 'missing_orthologs_counts.pdf'),
       width=6, height=4)

pg %>%
  filter(!grepl('^OG6.*', orthoMCL)) %>% 
  select(exp.name, Protein.IDs, Gene.names) %>% 
  write.table(file.path(out_folder, 'missing_orthologs_table.txt'),
              sep='\t', row.names=FALSE)

pg %>%
  select(exp.name, orthoMCL, id) %>% 
  mutate(orthoMCL=strsplit(orthoMCL, ';')) %>%
  unnest(orthoMCL) %>% 
  left_join(ortho_gene_names, by=c('orthoMCL'='ortho_group')) ->
  pg_gene_ids


pg %>%
  select(exp.name, Protein.IDs, id) %>% 
  mutate(Protein.IDs=strsplit(Protein.IDs, ';')) %>%
  unnest(Protein.IDs) %>%
  distinct_all() %>% 
  left_join(gene_names %>% select(Protein.IDs, Gene.names) %>% distinct_all(), 
            by='Protein.IDs') %>%   
  group_by(exp.name, id) %>% 
  summarise(specific.Gene.names=paste0(na.omit(unique(Gene.names)), collapse=';')) ->
  pg_gene_names

pg_gene_ids %>% 
  group_by(exp.name, id) %>%
  slice(1) %>%
  select(-orthoMCL) %>% 
  left_join(pg_gene_names, by=c('exp.name','id')) %>% 
  right_join(pg,
             by=c("exp.name","id")) %>% 
  mutate(plot_label=case_when(
    is.na(Gene.names_long) ~ sub('([^;]+);.*', '\\1;...', Protein.IDs),
    Gene.names_long == '' ~ sub('([^;]+);.*', '\\1;...', Protein.IDs), 
    TRUE ~ Gene.names_long)) %>% 
  ungroup() ->
  pg

pg %>% 
  select(orthoMCL, ortho_keywords, exp.name, Gene.names_long, Protein.IDs, specific.Gene.names) %>% 
  # select(exp.name, id, orthoMCL, ortho_keywords) %>%
  pivot_wider(names_from=exp.name, 
              values_from=c(Protein.IDs), 
              values_fn=~paste(.x, collapse=',')) ->
  ortho_group_by_species

write.table(ortho_group_by_species,
            file.path(out_folder, 'ortho_groups_by_species_full.txt'),
            sep='\t', row.names=FALSE)  



pg %>%
  group_by(exp.name) %>%
  summarise(n=n()) %>%
  ggplot(aes(exp.name, n)) +
  geom_col() + 
  coord_flip() +
  geom_text(aes(label=n), hjust=1.2, color='white') +
  labs(title='protein group count after filtering (2peps, 1unique)')
ggsave(file.path(out_folder,'proteingroup_counts.pdf'), width=11.69, height=6)

pg %>%
  select(exp.name, pulldown.date, Protein.IDs, id, orthoMCL, Sequence.coverage....,
         plot_label, specific.Gene.names, 
         starts_with('Identification.type'),
         starts_with('Intensity.'),
         starts_with('LFQ.intensity.'),
         starts_with('MS.MS.count.')) %>% 
  pivot_longer(c(starts_with('Intensity.'),
                 starts_with('LFQ.intensity.'),
                 starts_with('MS.MS.count.')), 
               names_to='col_names') %>%
  extract(col_names, 
          c('value_type', 'pulldown', 'replicate'),
          '(.*)\\.([^\\.]+)_(\\d)$') %>%
  pivot_wider(names_from=value_type, values_from=value) %>% 
  pivot_longer(starts_with('Identification.type'), 
               names_to='col_names') %>%
  extract(col_names, 
          c('value_type', 'pulldown2', 'replicate2'),
          '(.*)\\.([^\\.]+)_(\\d)$') %>%
  pivot_wider(names_from=value_type, values_from=value) %>%
  filter(replicate == replicate2, pulldown == pulldown2) %>%
  select(-replicate2, -pulldown2) %>%
  mutate(log2.Intensity=ifelse(Intensity == 0, NA, log2(Intensity)),
         log2.LFQ.intensity=ifelse(LFQ.intensity == 0, NA, log2(LFQ.intensity)),
         label=if_else(orthoMCL == '', Protein.IDs, orthoMCL)) %>% 
  ungroup() ->
  pg_long


pg_imputation_settings <- 
  pg_long %>%
  select(exp.name, pulldown, replicate, starts_with('log2')) %>%
  group_by(exp.name, pulldown, replicate) %>%
  summarise(n.imputed=sum(is.na(log2.LFQ.intensity)),
            n.Intensity=sum(!is.na(log2.Intensity)),
            min.log2.Intensity=min(log2.Intensity, na.rm=TRUE),
            max.log2.Intensity=max(log2.Intensity, na.rm=TRUE),
            mean.log2.Intensity=mean(log2.Intensity, na.rm=TRUE),
            median.log2.Intensity=median(log2.Intensity, na.rm=TRUE),
            sd.log2.Intensity=sd(log2.Intensity, na.rm=TRUE),
            low.log2.Intensity=quantile(log2.Intensity, low_percentile, na.rm=TRUE),
            high.log2.Intensity=quantile(log2.Intensity, high_percentile, na.rm=TRUE),
            n.LFQ.Intensity=sum(!is.na(log2.LFQ.intensity)),
            min.log2.LFQ.intensity=min(log2.LFQ.intensity, na.rm=TRUE),
            max.log2.LFQ.intensity=max(log2.LFQ.intensity, na.rm=TRUE),
            mean.log2.LFQ.intensity=mean(log2.LFQ.intensity, na.rm=TRUE),
            median.log2.LFQ.intensity=median(log2.LFQ.intensity, na.rm=TRUE),
            sd.log2.LFQ.intensity=sd(log2.LFQ.intensity, na.rm=TRUE),
            low.log2.LFQ.intensity=quantile(log2.LFQ.intensity, low_percentile, na.rm=TRUE),
            high.log2.LFQ.intensity=quantile(log2.LFQ.intensity, high_percentile, na.rm=TRUE))

pg_imputation_settings %>% 
  select(exp.name, pulldown, replicate, matches('log2.LFQ')) %>% 
  pivot_longer(matches('log2.LFQ'), names_to='type') %>% 
  mutate(type=sub('\\.log2.LFQ.intensity','', type)) %>% 
  ggplot(aes(exp.name, value, color=pulldown, shape=replicate)) + 
  geom_point(position=position_jitter(width=.3,height=0), alpha=.3) + 
  facet_wrap( ~ type, scale='free_y') + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title='comparing imputation settings')
ggsave(file.path(out_folder,'imputation_settings.pdf'), width=11.69, height=8.27)

set.seed(random_seed)

pg_long_imp <- 
  pg_long %>% 
  left_join(pg_imputation_settings,
            by=c('exp.name', 'pulldown', 'replicate')) %>% 
  mutate(imp.log2.LFQ.intensity=
           ifelse(is.na(log2.LFQ.intensity),
                  unlist(map2(low.log2.LFQ.intensity, high.log2.LFQ.intensity, ~ rbeta(1, 2, 2) * (.y - .x) + .x)),
                  log2.LFQ.intensity),
         imp.type=case_when(
           is.na(Intensity) ~ 'below detection',
           !is.na(Intensity) & is.na(log2.LFQ.intensity) ~ 'imputed',
           TRUE ~ 'quantified'
         )) %>%
  select(-names(pg_imputation_settings)[-1:-3])





pg_long_imp %>%
  group_by(exp.name, Protein.IDs, id, specific.Gene.names, label, plot_label, pulldown) %>%
  mutate(value_count=sum(!is.na(log2.LFQ.intensity))) %>%
  ungroup() %>%
  group_by(exp.name, Protein.IDs, label, plot_label, id, specific.Gene.names) %>%
  mutate(noLesion=rep(imp.log2.LFQ.intensity[pulldown == 'noLesion'], 4),
         value_count_noLesion=rep(value_count[pulldown == 'noLesion'], 4)) %>%
  ungroup() %>%
  filter(pulldown != 'noLesion') %>% 
  group_by(exp.name, Protein.IDs, label, plot_label, pulldown, id, specific.Gene.names) %>%
  summarise(p.value=t.test(imp.log2.LFQ.intensity, noLesion)$p.value,
            foldchange=mean(imp.log2.LFQ.intensity, na.rm=TRUE) - 
              mean(noLesion, na.rm=TRUE),
            value_count=unique(value_count),
            value_count_noLesion=unique(value_count_noLesion)) ->
  pg_avg


pg_long_imp %>% 
  ggplot(aes(pulldown, log10(LFQ.intensity), color=replicate)) +
  geom_boxplot() +
  facet_wrap(~exp.name)
ggsave(file.path(out_folder, 'intensity_boxplot.pdf'),
       width=11.69, height=8.27)



pg_long_imp %>% 
  select(exp.name, id, pulldown, replicate, log2.LFQ.intensity, imp.log2.LFQ.intensity) %>% 
  mutate(imputed=is.na(log2.LFQ.intensity)) %>% 
  filter(exp.name == 'A. thaliana') %>% 
  ggplot(aes(imp.log2.LFQ.intensity, fill=imputed, color=imputed)) + 
  labs(title='A. thaliana',
       subtitle='Distribution of imputed LFQ values') +
  geom_density(alpha=.4) +
  facet_grid(replicate ~ pulldown) +
  geom_rug() +
  scale_x_continuous(expand = c(0.01, 0)) +
  theme_minimal()
ggsave(file.path(out_folder, 'imputation_example.pdf'),
       width=11.69, height=8.27)



pg_avg %>%
  ungroup() %>%
  mutate(regulated=foldchange > 1 & p.value < 0.05,
    enriched=regulated & foldchange > 0
  ) ->
  pg_avg_enriched

pg_avg_enriched %>%
  left_join(pg_avg_enriched %>% 
              filter(enriched) %>%
              distinct(exp.name, label) %>%
              mutate(enriched_in_any=TRUE) %>% 
              rename(enriched_exp.name=exp.name)) %>%
  mutate(enriched_in_any=replace_na(enriched_in_any, FALSE)) ->
  pg_avg_enriched2
 

pg_avg_enriched2 %>%
  filter(enriched) %>% 
  group_by(label) %>% 
  distinct(exp.name, pulldown) ->
  jnk1 
jnk1 %>% 
  distinct(exp.name) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) ->
  jnk2
jnk1 %>% 
  distinct(pulldown) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) ->
  jnk3

jnk1 %>% 
  left_join(jnk1 %>% 
              distinct(exp.name) %>% 
              summarise(enriched_species=n())) %>% 
  left_join(jnk1 %>% 
              distinct(pulldown) %>% 
              summarise(enriched_pulldown=n())) %>% 
  right_join(pg_avg_enriched) %>% 
  rowwise() %>% 
  mutate(plot_label_long=paste(
    na.omit(c(plot_label, enriched_species, enriched_pulldown)), 
    collapse=',')) ->
  pg_avg_enriched


pg_avg_enriched %>% 
  filter(enriched) %>% 
  select(plot_label_long, exp.name, Protein.IDs) %>% 
  pivot_wider(names_from=exp.name, 
              values_from=c(Protein.IDs), 
              values_fn=~paste(.x, collapse=',')) ->
  ortho_group_by_species_enriched

write.table(ortho_group_by_species_enriched,
            file.path(out_folder,'ortho_groups_by_species_detected_enriched.txt'),
            sep='\t', row.names=FALSE)  

pg_avg_enriched %>% 
  filter(enriched) %>% 
  select(exp.name, label, plot_label, Protein.IDs) %>% 
  distinct() %>% 
  mutate(label=strsplit(label, ';')) %>% 
  unnest(label) %>% 
  left_join(orthoMCL %>% 
              select(ortho_group, ortho_Protein.IDs), 
            by=c('label'='ortho_group')) %>% 
  select(-exp.name, -Protein.IDs) %>% 
  extract(ortho_Protein.IDs, c('species','Protein.IDs'), '(.*)\\|(.*)') %>% 
  filter(!is.na(species)) %>% 
  pivot_wider(names_from=species, 
              values_from=Protein.IDs, 
              values_fn=~paste(.x, collapse=',')) ->
  ortho_group_by_species_full_enriched
write.table(ortho_group_by_species_full_enriched,
            file.path(out_folder,'ortho_groups_by_species_full_enriched.txt'),
            sep='\t', row.names=FALSE) 


pg_avg_enriched %>% 
  select(exp.name, label, plot_label, Protein.IDs) %>% 
  distinct() %>% 
  mutate(label=strsplit(label, ';')) %>% 
  unnest(label) %>% 
  left_join(orthoMCL %>% 
              select(ortho_group, ortho_Protein.IDs), 
            by=c('label'='ortho_group')) %>% 
  select(-exp.name, -Protein.IDs) %>% 
  extract(ortho_Protein.IDs, c('species','Protein.IDs'), '(.*)\\|(.*)') %>% 
  filter(!is.na(species)) %>% 
  pivot_wider(names_from=species, 
              values_from=Protein.IDs, 
              values_fn=~paste(.x, collapse=',')) ->
  ortho_group_by_species_full
write.table(ortho_group_by_species_full,
            file.path(out_folder,'ortho_groups_by_species_full.txt'),
            sep='\t', row.names=FALSE) 



pg_avg_enriched %>%
  filter(enriched_species == 1) %>% 
  select(exp.name, plot_label) %>% 
  distinct() %>% 
  inner_join(pg_avg_enriched) %>% 
  ggplot(aes(exp.name, plot_label, fill=enriched)) +
  geom_tile() +
  ggtitle('ortholog enriched in only one species') +
  facet_wrap(~pulldown, nrow=1) +
  scale_fill_manual(values=c(`TRUE`='#7f0000', `FALSE`='#AAAAAA')) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out_folder, 'single_enriched_ortholog.pdf'),
       width=8.27, height=23.38)



pg_avg_enriched %>%
  filter(
    enriched_pulldown >= 2 | 
      enriched_species >= 2, 
    enriched) %>% 
  select(exp.name, plot_label, Protein.IDs) %>% 
  distinct() %>% 
  inner_join(pg_avg_enriched) %>%
  group_by(label, plot_label, exp.name, pulldown) %>% 
  summarise_at(vars(enriched), any) %>% 
  pivot_wider(names_from=pulldown, values_from=enriched) ->
  pg_clustering


pg_clustering %>% 
  mutate(value=as.numeric(`8oxoG`)) %>% 
  select(-abasic, -RNAbase, -`8oxoG`) %>% 
  pivot_wider(names_from=exp.name, values_fill=-1) ->
  cluster_8oxoG

pg_clustering %>% 
  mutate(value=as.numeric(`RNAbase`)) %>% 
  select(-abasic, -RNAbase, -`8oxoG`) %>% 
  pivot_wider(names_from=exp.name, values_fill=-1) ->
  cluster_RNAbase

pg_clustering %>% 
  mutate(value=as.numeric(`abasic`)) %>% 
  select(-abasic, -RNAbase, -`8oxoG`) %>% 
  pivot_wider(names_from=exp.name, values_fill=-1) ->
  cluster_abasic

gplots::heatmap.2(
  cluster_8oxoG[3:14] %>% as.matrix(),
  Colv=TRUE, 
  labRow=unlist(cluster_8oxoG[,2])
)

gplots::heatmap.2(
  cluster_RNAbase[3:14] %>% as.matrix(),
  Colv=TRUE, 
  labRow=unlist(cluster_RNAbase[,2])
)

gplots::heatmap.2(
  cluster_abasic[3:14] %>% as.matrix(),
  Colv=TRUE, 
  labRow=unlist(cluster_abasic[,2])
)







pg_avg_enriched %>%
  filter(
    enriched_pulldown >= 2 |
      enriched_species >= 2,
    enriched) %>%
  select(exp.name, plot_label, Protein.IDs) %>%
  distinct() %>%
  inner_join(pg_avg_enriched) %>%
  group_by(label, plot_label, exp.name, pulldown) %>%
  summarise_at(vars(enriched), any) %>%
  # arrange(enriched) %>%
  ggplot(aes(exp.name, plot_label, fill=enriched)) +
  geom_tile() +
  ggtitle('ortholog enriched in minimum two species') +
  facet_wrap(~pulldown, nrow=1) +
  scale_fill_manual(values=c(`TRUE`='#7f0000', `FALSE`='#AAAAAA')) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out_folder, 'min2_enriched_ortholog.pdf'),
       width=8.3, height=11.7)

pg_avg_enriched %>%
  filter(
    enriched_pulldown >= 2 | 
      enriched_species >= 2, 
    enriched) %>% 
  select(exp.name, plot_label, Protein.IDs) %>% 
  distinct() %>% 
  inner_join(pg_avg_enriched) %>%
  group_by(label, plot_label, exp.name, pulldown) %>% 
  summarise_at(vars(enriched), any) %>%
  mutate(enriched=as.numeric(enriched)) %>% 
  pivot_wider(names_from=pulldown, values_from=enriched) %>% 
  pivot_wider(names_from=exp.name, 
              values_from=c(`8oxoG`, abasic, RNAbase),
              values_fill=-1) %>% 
  ungroup()->
  jnk
jnk %>% select(-label, -plot_label) %>% as.matrix() ->
  jnk2
rownames(jnk2) <- jnk$plot_label
gplots::heatmap.2(
  jnk2,
  Colv=FALSE, 
  dendrogram='row',
  margins=c(10,10)
)
jnk3 <- hclust(dist(jnk2))
jnk$plot_label[jnk3$order]

numeric_matrix <- jnk2
hclust_object <- jnk3
save(numeric_matrix, hclust_object, file='cluster_data.RData')
pdf(file.path(out_folder,'dendrogram.pdf'), 
    width=11.69, height=8.27) 
plot(hclust_object)
dev.off()

pg_avg_enriched %>%
  filter(
    # plot_label %in% c('zmay_100272929','zmay_100191897'),
    # plot_label %in% c('cele_cpl-1', 'zmay_100272929','zmay_100191897'),
    enriched_pulldown >= 2 |
      enriched_species >= 2,
    enriched) %>%
  select(exp.name, plot_label, Protein.IDs) %>%
  distinct() %>%
  inner_join(pg_avg_enriched) %>%
  group_by(label, plot_label, exp.name, pulldown) %>%
  summarise_at(vars(enriched), any) %>%
  mutate(plot_label=factor(plot_label, levels=jnk$plot_label[rev(jnk3$order)])) %>% 
  # arrange(enriched) %>%
  ggplot(aes(exp.name, plot_label, fill=enriched)) +
  geom_tile() +
  labs(title='ortholog enriched in minimum two species',
       subtitle='clustered by hclust on numerical matrix c(-1,0,1)') +
  facet_wrap(~pulldown, nrow=1) +
  scale_fill_manual(values=c(`TRUE`='#7f0000', `FALSE`='#AAAAAA')) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out_folder, 'cluster_min2_enriched_ortholog.pdf'),
              width=8.3, height=11.7)



pg_avg_enriched %>%
  filter(
    # plot_label %in% c('zmay_100272929','zmay_100191897'),
    # plot_label %in% c('cele_cpl-1', 'zmay_100272929','zmay_100191897'),
    enriched_pulldown >= 2 |
      enriched_species >= 2,
    enriched) %>%
  select(exp.name, plot_label, Protein.IDs) %>%
  distinct() %>%
  inner_join(pg_avg_enriched) %>%
  group_by(label, plot_label, exp.name, pulldown) %>%
  summarise(foldchange=max(foldchange)) %>%
  mutate(plot_label=factor(plot_label, levels=jnk$plot_label[rev(jnk3$order)])) %>% 
  # arrange(enriched) %>%
  ggplot(aes(exp.name, plot_label, fill=foldchange)) +
  geom_tile() +
  labs(title='ortholog enriched in minimum two species',
       subtitle='clustered by hclust on numerical matrix c(-1,0,1)') +
  facet_wrap(~pulldown, nrow=1) +
  scale_fill_gradient('',
                      low='#ffeda0', high='#800026', 
                      limits=c(0,1), 
                      breaks=c(0,log2(1.25),log2(1.5),log2(1.75),1),
                      labels=c('not enriched','1.25-fold','1.50-fold','1.75-fold','enriched'),
                      oob = scales::squish) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out_folder, 'cluster_min2_gradient_ortholog.pdf'),
       width=8.3, height=11.7)



pg_avg_enriched %>%
  filter(
    # plot_label %in% c('zmay_100272929','zmay_100191897'),
    # plot_label %in% c('cele_cpl-1', 'zmay_100272929','zmay_100191897'),
    enriched_pulldown >= 2 &
      enriched_species == 1,
    enriched) %>%
  select(exp.name, plot_label, Protein.IDs) %>%
  distinct() %>%
  inner_join(pg_avg_enriched) %>%
  group_by(label, plot_label, exp.name, pulldown) %>%
  summarise(foldchange=max(foldchange),
            specific.Gene.names=ifelse(unique(specific.Gene.names) == '',
                                       plot_label, unique(specific.Gene.names))) %>%
  mutate(plot_label=factor(plot_label, levels=jnk$plot_label[rev(jnk3$order)])) %>% 
  arrange(plot_label) ->
  jnk_heatmap
  # arrange(enriched) %>%
jnk_heatmap %>% 
  ggplot(aes(exp.name, plot_label, fill=foldchange)) +
  geom_tile() +
  labs(title='ortholog enriched in one organism, multiple legions',
       subtitle='clustered by hclust on numerical matrix c(-1,0,1)') +
  facet_wrap(~pulldown, nrow=1) +
  scale_fill_gradient('',
                      low='#ffeda0', high='#800026', 
                      limits=c(0,1), 
                      breaks=c(0,log2(1.25),log2(1.5),log2(1.75),1),
                      labels=c('not enriched','1.25-fold','1.50-fold','1.75-fold','enriched'),
                      oob = scales::squish) +
  # scale_y_discrete(label=unique(jnk_heatmap$specific.Gene.names)) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out_folder, 'cluster_min2_gradient_ortholog_1species_multipleLegions.pdf'),
       width=8.3, height=11.7)
pg_avg_enriched %>%
  filter(
    # plot_label %in% c('zmay_100272929','zmay_100191897'),
    # plot_label %in% c('cele_cpl-1', 'zmay_100272929','zmay_100191897'),
    enriched_pulldown >= 1 &
      enriched_species >= 2,
    enriched) %>%
  select(exp.name, plot_label, Protein.IDs) %>%
  distinct() %>%
  inner_join(pg_avg_enriched) %>%
  group_by(label, plot_label, exp.name, pulldown) %>%
  summarise(foldchange=max(foldchange)) %>%
  mutate(plot_label=factor(plot_label, levels=jnk$plot_label[rev(jnk3$order)])) %>% 
  # arrange(enriched) %>%
  ggplot(aes(exp.name, plot_label, fill=foldchange)) +
  geom_tile() +
  labs(title='ortholog enriched in minimum two species and min one legion',
       subtitle='clustered by hclust on numerical matrix c(-1,0,1)') +
  facet_wrap(~pulldown, nrow=1) +
  scale_fill_gradient('',
                      low='#ffeda0', high='#800026', 
                      limits=c(0,1), 
                      breaks=c(0,log2(1.25),log2(1.5),log2(1.75),1),
                      labels=c('not enriched','1.25-fold','1.50-fold','1.75-fold','enriched'),
                      oob = scales::squish) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out_folder, 'cluster_min2_gradient_ortholog_1Legions_2species.pdf'),
       width=8.3, height=11.7)


my_special_dendrogram <- 
  c("hsap_APEX1", "scer_APN1", "hsap_CRY1", "hsap_RFC3", "hsap_RFC4", "hsap_RFC2", 
    "hsap_RFC1", "scer_ASG1", "hsap_MYO15A", "hsap_MYBBP1A", "scer_RSC58", 
    "scer_PDR3", "scer_SWI6", "hsap_WDR76", "hsap_INO80", "hsap_PBRM1", 
    "hsap_SMARCD1", "hsap_SMARCA2", "hsap_TOP2B", "hsap_MUTYH", "hsap_MPG", 
    "bsub_ydaT", "bsub_yfjM", "bsub_ydeI", "ecol_dinG", "hsap_TOP3B", 
    "ecol_mutM", "hsap_KCTD9", "bsub_yhaZ", "hsap_OGG1", "hsap_POLN", 
    "cele_CELE_F07A5.2", "cele_CELE_T01E8.8", "cele_col-64", 
    "tetr_TTHERM_00681750", "cele_CELE_Y75B7B.2", "SPAC3H8.08c", 
    "A0A1D6LV91;...", "zmay_100191897", "A0A1D6NSE6", "zmay_100192956", 
    "zmay_GBP20", "zmay_100272929", "cele_hmg-12", "atha_CRYD", "hsap_TREH", 
    "atha_MOC1", "atha_PHR1", "spom_hpz2", "hsap_RBM34", "hsap_RAE1", 
    "hsap_PNKP", "hsap_DNAJC13", "hsap_APTX", "hsap_LIG3", "hsap_PARP2", 
    "hsap_POLB", "hsap_XRCC1", "hsap_NIP7", "hsap_CENPV", "hsap_CHD1", 
    "hsap_KAT6A", "hsap_AHCTF1", "hsap_PCGF1", "hsap_ZNF512B", "hsap_NOC3L", 
    "hsap_MYL9")


pg_avg_enriched %>%
  filter(
    # plot_label %in% c('zmay_100272929','zmay_100191897'),
    # plot_label %in% c('cele_cpl-1', 'zmay_100272929','zmay_100191897'),
    enriched_pulldown >= 1 &
      enriched_species >= 2,
    enriched) %>%
  select(exp.name, plot_label, Protein.IDs) %>%
  distinct() %>%
  inner_join(pg_avg_enriched) %>%
  group_by(label, plot_label, exp.name, pulldown) %>%
  summarise(foldchange=max(foldchange)) %>%
  mutate(plot_label=factor(plot_label, levels=rev(my_special_dendrogram))) %>% 
  # arrange(enriched) %>%
  ggplot(aes(exp.name, plot_label, fill=foldchange)) +
  geom_tile() +
  labs(title='ortholog enriched in minimum two species and min one legion',
       subtitle='clustered by hclust on numerical matrix c(-1,0,1)') +
  facet_wrap(~pulldown, nrow=1) +
  scale_fill_gradient('',
                      low='#ffeda0', high='#800026', 
                      limits=c(0,1), 
                      breaks=c(0,log2(1.25),log2(1.5),log2(1.75),1),
                      labels=c('not enriched','1.25-fold','1.50-fold','1.75-fold','enriched'),
                      oob = scales::squish) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out_folder, 'cluster_min2_gradient_ortholog_1Legions_2species_manualOrder.pdf'),
       width=8.3, height=11.7)





set.seed(1)
pg_avg_enriched %>% 
  ungroup() %>%
  distinct(foldchange, p.value, enriched_species, enriched_pulldown, enriched,
           exp.name, pulldown, plot_label_long) %>% 
  # slice_sample(n=10000) %>%
  mutate(plot_label_long=if_else(
    enriched, 
    sub('([^;]+);.*', '\\1;...', plot_label_long), 
    '')) %>% 
  ggplot(aes(foldchange, -log10(p.value), 
             shape=enriched_species>1, size=enriched_pulldown, fill=enriched, 
             color=enriched_species>0&!enriched)) +
  geom_point(data=. %>% filter(is.na(enriched_species)), size=.7, shape=1, alpha=.5) +
  geom_hline(yintercept=-log10(0.05),linetype='dashed') +
  geom_vline(xintercept=1, linetype='dashed') +
  geom_point(data=. %>% filter(!is.na(enriched_species))) +
  facet_wrap(exp.name ~ pulldown, scales='free') +
  scale_shape_manual('enriched in more', values=c('FALSE'=21, 'TRUE'=22)) +
  scale_fill_manual(values=c('FALSE'='#99999955', 'TRUE'='#FF0000FF')) +
  scale_size_continuous(range=c(1,2)) +
  scale_color_manual(values=c('FALSE'='#99999955', 'TRUE'='#FF000088')) +
  geom_text_repel(aes(label=plot_label_long), max.overlaps=100, size=2, color='#FF0000FF') +
  scale_x_continuous(limits=c(-3,max(pg_avg_enriched$foldchange))) +
  scale_y_continuous(limits=c(0,max(-log10(pg_avg_enriched$p.value)))) +
  theme_light(7)
ggsave(file.path(out_folder,'volcano_all.pdf'), 
       width=2*11.69, height=2*8.27) 

set.seed(1)
pg_avg_enriched %>% 
  ungroup() %>%
  distinct(foldchange, p.value, enriched_species, enriched_pulldown, enriched,
           exp.name, pulldown, plot_label_long) %>%
  mutate(plot_label_long=if_else(
    enriched, 
    sub('([^;]+);.*', '\\1;...', plot_label_long), 
    '')) %>% 
  ggplot(aes(foldchange, -log10(p.value), 
             shape=enriched_species>1, size=enriched_pulldown, fill=enriched, 
             color=enriched_species>0&!enriched)) +
  geom_hline(yintercept=-log10(0.05),linetype='dashed') +
  geom_vline(xintercept=1, linetype='dashed') +
  geom_point(data=. %>% filter(!is.na(enriched_species))) +
  facet_wrap(exp.name ~ pulldown, scales='free') +
  scale_shape_manual('enriched in more', values=c('FALSE'=21, 'TRUE'=22)) +
  scale_fill_manual(values=c('FALSE'='#99999955', 'TRUE'='#FF0000FF')) +
  scale_size_continuous(range=c(1,2)) +
  scale_color_manual(values=c('FALSE'='#99999955', 'TRUE'='#FF000088')) +
  geom_text_repel(aes(label=plot_label_long), max.overlaps=100, size=2, color='#FF0000FF') +
  scale_x_continuous(limits=c(-1,max(pg_avg_enriched$foldchange))) +
  scale_y_continuous(limits=c(0,max(-log10(pg_avg_enriched$p.value)))) +
  theme_light(7)
ggsave(file.path(out_folder,'volcano_enriched_any.pdf'), 
       width=2*11.69, height=2*8.27) 


set.seed(1)
pg_avg_enriched %>% 
  ungroup() %>%
  distinct(foldchange, p.value, enriched_species, id, enriched_pulldown, enriched,
           exp.name, pulldown, plot_label) %>% 
  left_join(pg_avg_enriched %>% 
              ungroup() %>% 
              filter(enriched) %>% 
              select(exp.name, pulldown, id, plot_label) %>% 
              group_by(plot_label, pulldown) %>% 
              mutate(enriched_advanced=if_else(length(unique(exp.name)) >= 2, 
                                               'same legion; >=2 species', 
                                               'single species')) %>% 
              ungroup() %>%  
              select(-plot_label)) %>% 
  mutate(plot_label_long=if_else(
    enriched, 
    sub('([^;]+);.*', '\\1;...', plot_label), 
    ''),
    enriched_advanced=if_else(is.na(enriched_advanced), 'background', enriched_advanced)) %>% 
  ggplot(aes(foldchange, -log10(p.value), 
             # fill=enriched, 
             color=enriched_advanced)) +
  geom_point(data=. %>% filter(is.na(enriched_species)),alpha=.5) +
  geom_hline(yintercept=-log10(0.05),linetype='dashed') +
  geom_vline(xintercept=1, linetype='dashed') +
  geom_point(data=. %>% filter(!is.na(enriched_species))) +
  facet_wrap(exp.name ~ pulldown, scales='free') +
  scale_color_manual(values=c('background'='#99999955', 
                              'same legion; >=2 species'='#0000FF88', 
                              'single species'='#FF000088')) +
  geom_text_repel(aes(label=plot_label_long), max.overlaps=100, size=2) +
  scale_x_continuous(limits=c(-3,max(pg_avg_enriched$foldchange))) +
  scale_y_continuous(limits=c(0,max(-log10(pg_avg_enriched$p.value)))) +
  theme_light(7)
ggsave(file.path(out_folder,'volcano_new_color.pdf'), 
       width=2*11.69, height=2*8.27) 


set.seed(1)
pg_avg_enriched %>% 
  ungroup() %>%
  distinct(foldchange, p.value, enriched_species, id, enriched_pulldown, enriched,
           exp.name, pulldown, plot_label, specific.Gene.names, Protein.IDs) %>% 
  mutate(plot_label_long=if_else(
    enriched, 
    if_else(exp.name == 'Z. mays',
            sub('([^;]+);.*', '\\1', Protein.IDs),
            if_else(is.na(specific.Gene.names) | specific.Gene.names == '',
                    sub('([^;]+);.*', '\\1', Protein.IDs),
                    sub('([^;]+);.*', '\\1', specific.Gene.names))
    ), 
    ''),
    enriched_advanced=if_else(enriched, 'enriched','background')
    ) %>% 
  ggplot(aes(foldchange, -log10(p.value), 
             # fill=enriched,
             color=enriched_advanced)) +
  geom_point(data=. %>% filter(is.na(enriched_species)),alpha=.5) +
  geom_hline(yintercept=-log10(0.05),linetype='dashed') +
  geom_vline(xintercept=1, linetype='dashed') +
  geom_point(data=. %>% filter(!is.na(enriched_species))) +
  facet_wrap(exp.name ~ pulldown, scales='free') +
  scale_color_manual(values=c('background'='#99999955',
                              'enriched'='#0000FF88')) +
  geom_text_repel(aes(label=plot_label_long), max.overlaps=1000, 
                  size=1.5, color='black') +
  scale_x_continuous(limits=c(-3,max(pg_avg_enriched$foldchange))) +
  scale_y_continuous(limits=c(0,max(-log10(pg_avg_enriched$p.value)))) +
  theme_light(7) ->
  volcano_plot
ggsave(file.path(out_folder,'volcano_speciesGeneNames.pdf'), 
       volcano_plot,
       width=2*11.69, height=2*8.27) 

set.seed(1)
pg_avg_enriched %>% 
  ungroup() %>%
  distinct(foldchange, p.value, enriched_species, enriched_pulldown, enriched,
           exp.name, pulldown, plot_label_long, Protein.IDs, specific.Gene.names) %>%
  mutate(plot_label_long=if_else(
    enriched, 
    if_else(exp.name == 'Z. mays',
            sub('([^;]+);.*', '\\1', Protein.IDs),
            if_else(is.na(specific.Gene.names) | specific.Gene.names == '',
                    sub('([^;]+);.*', '\\1', Protein.IDs),
                    sub('([^;]+);.*', '\\1', specific.Gene.names))
    ), 
    # ''),
    as.character(NA)),
    enriched_advanced=if_else(enriched, 'enriched','background')) %>%
  ungroup() %>% 
  group_by(pulldown) %>% 
  do(p=ggplot(data=., 
              aes(foldchange, -log10(p.value), 
                  color=enriched_advanced)) +
       geom_point(data=. %>% filter(is.na(enriched_species)), size=.7) +
       geom_hline(yintercept=-log10(0.05),linetype='dashed') +
       geom_vline(xintercept=1, linetype='dashed') +
       geom_point(data=. %>% filter(!is.na(enriched_species))) +
       facet_wrap(~ exp.name, scales='free') +
       ggtitle(unique(.$pulldown)) +
       scale_color_manual(values=c('background'='#99999955',
                                   'enriched'='#0000FFAA')) +
       geom_text_repel(aes(label=plot_label_long), max.overlaps=100, 
                       size=2, color='black') +
       scale_x_continuous(limits=c(-3,max(pg_avg_enriched$foldchange))) +
       scale_y_continuous(limits=c(0,max(-log10(pg_avg_enriched$p.value)))) +
       theme_light(7)) ->
  plots
plots$p[[1]]
pdf(file.path(out_folder,'volcano_perDamage_speciesGeneNames.pdf'), 
    width=11.69, height=8.27) 
plots$p
dev.off()

pg_avg_enriched %>% 
  ungroup() %>%
  distinct(foldchange, p.value, enriched_species, enriched_pulldown, enriched,
           exp.name, pulldown, plot_label_long) %>%
  # slice_sample(n=10000) %>%
  mutate(plot_label_long=if_else(
    enriched, 
    sub('([^;]+);.*', '\\1;...', plot_label_long), 
    '')) %>%
  ungroup() %>% 
  group_by(pulldown) %>% 
  do(p=ggplot(data=., 
              aes(foldchange, -log10(p.value), 
                  shape=enriched_species>1, size=enriched_pulldown, fill=enriched, 
                  color=enriched_species>0&!enriched)) +
       geom_point(data=. %>% filter(is.na(enriched_species)), size=.7, shape=1, alpha=.5) +
       geom_hline(yintercept=-log10(0.05),linetype='dashed') +
       geom_vline(xintercept=1, linetype='dashed') +
       geom_point(data=. %>% filter(!is.na(enriched_species))) +
       facet_wrap(~ exp.name, scales='free') +
       scale_shape_manual('enriched in more', values=c('FALSE'=21, 'TRUE'=22)) +
       scale_fill_manual(values=c('FALSE'='#99999955', 'TRUE'='#FF0000FF')) +
       scale_size_continuous(range=c(1,2)) +
       ggtitle(unique(.$pulldown)) +
       scale_color_manual(values=c('FALSE'='#99999955', 'TRUE'='#FF000088')) +
       geom_text_repel(aes(label=plot_label_long), max.overlaps=100, size=2, color='#FF0000FF') +
       scale_x_continuous(limits=c(-3,max(pg_avg_enriched$foldchange))) +
       scale_y_continuous(limits=c(0,max(-log10(pg_avg_enriched$p.value)))) +
       theme_light(7)) ->
  plots

pdf(file.path(out_folder,'volcano_perDamage.pdf'), 
    width=11.69, height=8.27) 
plots$p
dev.off()



pg_avg_enriched %>% 
  ungroup() %>%
  distinct(foldchange, p.value, enriched_species, enriched_pulldown, enriched,
           exp.name, pulldown, plot_label_long) %>%
  mutate(plot_label_long=if_else(
    enriched, 
    sub('([^;]+);.*', '\\1;...', plot_label_long), 
    '')) %>%
  ungroup() %>% 
  group_by(exp.name) %>% 
  do(p=ggplot(data=., 
              aes(foldchange, -log10(p.value), 
                  shape=enriched_species>1, size=enriched_pulldown, fill=enriched, 
                  color=enriched_species>0&!enriched)) +
       geom_point(data=. %>% filter(is.na(enriched_species)), size=.7, shape=1, alpha=.5) +
       geom_hline(yintercept=-log10(0.05),linetype='dashed') +
       geom_vline(xintercept=1, linetype='dashed') +
       geom_point(data=. %>% filter(!is.na(enriched_species))) +
       facet_wrap(~ pulldown, scales='free') +
       scale_shape_manual('enriched in more', values=c('FALSE'=21, 'TRUE'=22)) +
       scale_fill_manual(values=c('FALSE'='#99999955', 'TRUE'='#FF0000FF')) +
       scale_size_continuous(range=c(1,2)) +
       ggtitle(unique(.$exp.name)) +
       scale_color_manual(values=c('FALSE'='#99999955', 'TRUE'='#FF000088')) +
       geom_text_repel(aes(label=plot_label_long), max.overlaps=100, size=2, color='#FF0000FF') +
       # scale_x_continuous(limits=range(pg_avg_enriched$foldchange)) +
       scale_x_continuous(limits=c(-3,max(pg_avg_enriched$foldchange))) +
       scale_y_continuous(limits=c(0,max(-log10(pg_avg_enriched$p.value)))) +
       theme_light(7)) ->
  plots

pdf(file.path(out_folder,'volcano_perSpecies.pdf'), 
    width=11.69, height=8.27) 
plots$p
dev.off()


rMQanalysis::write.table_imb(
  as.data.frame(pg_avg_enriched), 
  file.path(out_folder, 'pg_avg_enriched.txt'))


pg_avg_enriched %>% 
  filter(enriched) %>%
  distinct(exp.name, label, pulldown) %>% 
  mutate(values=1) %>% 
  pivot_wider(names_from=pulldown, values_from=values) %>% 
  select(-label) ->
  jnk
pdf(file.path(out_folder,'venn_perSpecies.pdf'), 
    width=11.69, height=8.27) 
# par(mfrow = c(4, 3))
for(exp_name in unique(jnk$exp.name)) {
  jnk %>% 
    filter(exp.name == exp_name) %>%
    ungroup() %>%
    select(-exp.name) %>% 
    mutate_all(replace_na, 0) %>% 
    eulerr::euler() ->
    jnk2
    print(plot(jnk2, main=exp_name, counts = T, fontface = 1, quantities = TRUE))
}
dev.off()



