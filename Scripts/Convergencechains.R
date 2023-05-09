# plot chain
chain_conv = function(burnreps, df_chain=mainmod$chain,sfix){
  df_chain= df_chain%>% as.data.frame()
  names(df_chain) =c("Beta1","Beta2","Beta3","omega")
  df_chain = df_chain%>% drop_na()
  df_chain <- df_chain[burnreps:nrow(df_chain),]
  df_chain$MCMC = rownames(df_chain)
  df_chain = df_chain %>% gather(p,value,1:4)
  df_chain$MCMC = gsub("MCMC ","",df_chain$MCMC) %>% as.numeric
  p1=ggplot(df_chain, aes(x=MCMC, y=value, group=p))+
    geom_line() +
    facet_wrap(~p,scales="free_y")
  return(p1)
}

chain_conv2 = function(burnreps, df_chain=mainmod2$chain,sfix){
  df_chain= df_chain%>% as.data.frame()
  names(df_chain) =c("Beta1","Beta2","Beta3","omega")
  df_chain = df_chain%>% drop_na()
  df_chain <- df_chain[burnreps:nrow(df_chain),]
  df_chain$MCMC = rownames(df_chain)
  df_chain = df_chain %>% gather(p,value,1:4)
  df_chain$MCMC = gsub("MCMC ","",df_chain$MCMC) %>% as.numeric
  p1=ggplot(df_chain, aes(x=MCMC, y=value, group=p))+
    geom_line() +
    facet_wrap(~p,scales="free_y")
  return(p1)
}