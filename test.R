

# rand_object<-RandCapGen(
#   n=150,
#   arms = c("A","B","C"),
#   ratio=c(1,2,3),
#   block_sizes = c(6,12),
#   strat_vars = list(
#     age=c(1,2,3),
#     sex=c(1,2),
#     az=c(1,0),
#     centre=c(700:705)
#   ),
#   strat_vars_prefix = c(
#     age="age",
#     sex="sex",
#     az="az",
#     centre="c"
#   ),
#   seed = 2
# )


rand_object<-RandCapGen(
  n=240,
  arms = c("A","B"),
  block_sizes = c(4,6),
  project_acronym = "ANTEIPA"
)

RandCapBalance(rand_object)
RandCapSettings(rand_object)
prod_list<-RandCapProd(rand_object,5)
RandCapBalance(prod_list)
RandCapSettings(prod_list)
RandCapTable(prod_list,save_for_REDCap = T)

colnames(prod_list$tables$simplified_dataset)<-c(
  "rando_bras", "rando_age", "rando_sex", "rando_az", "redcap_data_access_group"
)

prod_list$tables$full_dataset%>%View()

