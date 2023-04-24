
### All functions for exposure_ratio_analyses

# 1. each column gets multiplied by a row in sharing matrix 
clean_ls = function(l_mat){
  clean_l = l_mat %>% 
    dplyr::select(Host:value) %>% 
    arrange(Host, name) %>% 
    pivot_wider(names_from=name, values_from=value) %>% 
    dplyr::rename(sp1=Host)
  return(clean_l)
}


# 2. find shared parasites
calculate_shared_parasites = function(sharing, data_template){
  sharing2 = rbind(sharing, sharing, sharing, sharing, sharing)
  sharing3 = rbind(sharing2, sharing2, sharing2)
  sharing3 = rbind(sharing3, sharing3, sharing3)
  sharing3 = rbind(sharing3, sharing3)
  sharing3 = sharing3 %>% 
    arrange(sp1)
  
  # check same dimensions
  if(length(sharing3$sp1)!=length(data_template$focal_sp)){
    print("Dimensions not compatible")
    print(dim(sharing3))
    print(dim(data_template))
  }
  
  opts = data_template %>% 
    ungroup() %>% 
    dplyr::select(BUFFALO:ZEBRA)
  
  new_parasite_mat = opts*sharing3[,-1]
  
  new_parasite_mat2 = cbind(data_template[,c(1:4,11)], new_parasite_mat)
  new_parasite_mat2$new_parasites = rowSums(new_parasite_mat2[,-c(1:6)]) 
  return(new_parasite_mat2)
}


# 3. reformat names
clean_output = function(out, oprisk){
  new_oprisk_a = out %>% 
    mutate(focal_sp=recode(focal_sp,
                           Buffalo="BUFFALO",
                           Cattle = "CATTLE",
                           Elephant = "ELEPHANT",
                           Impala = "IMPALA",
                           Giraffe = "GIRAFFE",
                           Zebra="ZEBRA")) %>% 
    dplyr::rename("choice"=focal_sp) %>% 
    left_join(oprisk, multiple="all")
}


# 4. filter to relevant rows and activity type
format_join_to_exposure = function(join_df){
  new_df = join_df %>%
    mutate_at(vars(BUFFALO:ZEBRA, parasite_density), funs(.*Ind_Secs/100)) %>%  ### divide to get per m2 (approximate area of camera)
    filter(Treatment %in% c("WPC","CONT"), Activity=="dgraze") %>%
    pivot_longer(BUFFALO:ZEBRA, names_to="other_sp",values_to="exposure") %>%
    group_by(Period, Site, Treatment, choice, other_sp) %>% 
    summarize_at(vars(exposure), funs(mean(., na.rm=T))) %>% 
    mutate(other_sp = factor(other_sp,levels=c("BUFFALO",
                                               "CATTLE",
                                               "GIRAFFE",
                                               "IMPALA",
                                               "ZEBRA",
                                               "ELEPHANT"))) %>% 
    mutate(choice = factor(choice, levels=c("BUFFALO",
                                            "CATTLE",
                                            "GIRAFFE",
                                            "IMPALA",
                                            "ZEBRA",
                                            "ELEPHANT"))) %>% 
    mutate(Treatment = recode(Treatment, "CONT"="Matrix","WPC"="Water")) 
  return(new_df)
}



# 5. Fit GLMMs for each species separately
fit_species_specific_exposure_model = function(species, df){
  
  # exclude species with zero exposure as this greatly inflates zeros and results in poor fit
  to_exclude = df %>% 
    group_by(choice, other_sp) %>% 
    summarize_at(vars(exposure), funs(sum(.,na.rm=T))) %>% 
    filter(exposure==0) %>% 
    filter(choice==species)
  
  l_other_sp = levels(as.factor(df$other_sp))
  l_other_sp = l_other_sp[-which(l_other_sp == species)]
  
  b_mepd3 = df %>% 
    filter(choice==species) %>% 
    filter(is.na(exposure)==F) %>% 
    filter(!c(other_sp %in% to_exclude$other_sp)) %>% 
    mutate(other_sp = factor(other_sp, levels=c(species, l_other_sp)))
  
  
  if(length(levels(as.factor(as.character(b_mepd3$other_sp))))>1){
    bmod = glmmTMB(exposure~other_sp*Treatment+(1|Site)+(1|Period), data=b_mepd3, family="tweedie")
    testResiduals(simulateResiduals(bmod))
    options(scipen=100000)
    
    # posthoc comparisons and formatting
    b_p_c = pairs(emmeans(bmod, ~other_sp|Treatment), adjust="none", type="response") %>%
      as.data.frame() %>% 
      as.data.frame() %>% 
      filter(grepl(species, contrast)==T)
    b_p_c$p.adjusted=p.adjust(b_p_c$p.value, method="holm")
    b_p_c$choice=species
    b_p_c$other_sp = sub(species, "", b_p_c$contrast)
    b_p_c$other_sp = sub(" / ", "", b_p_c$other_sp)
    
    b_e_c = emmeans(bmod, ~other_sp|Treatment, type="response") %>%
      as.data.frame() %>%
      as.data.frame()
    b_e_c$choice=species
  }else{
    b_p_c = 0
    b_e_c = 0
  }
  
  ## excluded
  b_mepd4 = df %>% 
    filter(choice==species) %>% 
    filter(is.na(exposure)==F) %>% 
    filter(c(other_sp==choice))
  
  if(length(levels(as.factor(as.character(b_mepd3$other_sp))))==1){
    bmod2 = glmmTMB(exposure~Treatment+(1|Site)+(1|Period), data=b_mepd4, family="tweedie")
    b_e_c2 = emmeans(bmod2, ~Treatment, type="response") %>%
      as.data.frame() %>%
      as.data.frame()
    b_e_c2$choice=species
    b_e_c2$other_sp=species
    b_e_c2 = b_e_c2 %>% 
      relocate(other_sp)
  }else{
    b_e_c2 = 0
  }
  
  
  return(list(pairs = b_p_c, means = b_e_c, extra = b_e_c2))
}



# 6. Create plot for each species' results
format_and_plot_sp_model_results = function(all_p_c, all_e_c, exp_plot_dat){
  # format for plot
  all_p_c = all_p_c %>% 
    mutate(sig = ifelse(p.adjusted <0.001, "***",
                        ifelse(p.adjusted <0.01, "**",
                               ifelse(p.adjusted <0.05, "*",
                                      ifelse(p.adjusted <0.1, "."," ")))))
  all_p_c = all_p_c %>% 
    dplyr::select(Treatment, p.adjusted, choice, other_sp, sig)
  add = data.frame(Treatment = rep(c("Matrix","Water"),4),
                   p.adjusted=NA,
                   choice=c(rep("BUFFALO",2),rep("CATTLE",2),rep("GIRAFFE",2),rep("IMPALA",2)),
                   other_sp=c(rep("BUFFALO",2),rep("CATTLE",2),rep("GIRAFFE",2),rep("IMPALA",2)),
                   sig=" ")
  
  all_p_c = rbind(all_p_c,add)
  
  maxes = exp_plot_dat %>% 
    group_by(Treatment, choice, other_sp) %>% 
    summarize_at(vars(exposure), funs(max(., na.rm=T)))
  all_p_c = left_join(all_p_c, maxes)
  all_p_c = exp_plot_dat %>% 
    group_by(choice, other_sp, Treatment) %>% 
    dplyr::summarize(n()) %>% 
    left_join(all_p_c) %>% 
    mutate_at(vars(sig), funs(ifelse(is.na(.), " ", .)))
  all_p_c = all_p_c %>% 
    mutate_at(vars(choice, other_sp), funs(str_to_title(.))) %>% 
    mutate(choice=factor(choice, levels=c("Cattle","Buffalo","Giraffe","Impala","Zebra","Elephant"))) %>% 
    mutate(other_sp=factor(other_sp, levels=c("Cattle","Buffalo","Giraffe","Impala","Zebra","Elephant")))
  
  # points
  exp_plot_dat2 = exp_plot_dat %>% 
    group_by(Site, Treatment, choice, other_sp) %>% 
    summarize_at(vars(exposure), funs(mean(., na.rm=T))) %>% 
    mutate_at(vars(choice, other_sp), funs(str_to_title(.))) %>% 
    mutate(choice=factor(choice, levels=c("Cattle","Buffalo","Giraffe","Impala","Zebra","Elephant"))) %>% 
    mutate(other_sp=factor(other_sp, levels=c("Cattle","Buffalo","Giraffe","Impala","Zebra","Elephant")))
  
 # errorbars
  add_sp = all_e_cm %>% 
    group_by(choice) %>% 
    dplyr::summarize(n=n()) %>% 
    filter(n!=8)
  
  if(length(add_sp$choice>0)){
    add2 = data.frame(other_sp = rep(add_sp$choice, 8),
                      Treatment = rep(rep(c("Matrix","Water"), each=length(add_sp$choice)), 4),
                      response=NA,
                      SE=NA,
                      df=NA,
                      lower.CL=NA,
                      upper.CL=NA,
                      choice = rep(c("Buffalo","Cattle","Giraffe","Impala"), each = length(add_sp$choice)*2))
    
    add3 = data.frame(other_sp = rep( rep(c("Buffalo","Cattle","Giraffe","Impala"), each = 2,), length(add_sp$choice)),
                      Treatment = rep(c("Matrix","Water"), each = length(add_sp$choice)*2),
                      response=NA,
                      SE=NA,
                      df=NA,
                      lower.CL=NA,
                      upper.CL=NA,
                      choice = c(rep(add_sp$choice, each=8)))
    
    all_e_c = rbind(all_e_c, add2, add3)
  }
  
  all_e_c = all_e_c %>% 
    mutate_at(vars(choice, other_sp), funs(str_to_title(.))) %>% 
    mutate(choice=factor(choice, levels=c("Cattle","Buffalo","Giraffe","Impala","Zebra","Elephant"))) %>% 
    mutate(other_sp=factor(other_sp, levels=c("Cattle","Buffalo","Giraffe","Impala","Zebra","Elephant"))) 
  
  
  exp_plot = {ggplot(exp_plot_dat2, aes(x=Treatment, y=(exposure+1)))+
     
      geom_point(aes(col=other_sp), position=position_dodge(width=0.5), size=1, alpha=0.5, shape = 15)+
      geom_point(data=all_e_c, aes(x=Treatment, y=response+1, col=other_sp), position=position_dodge(width=0.5), size=2)+
      geom_errorbar(data=all_e_c, aes(x=Treatment, y=response+1, ymin=lower.CL+1, ymax=upper.CL+1, col=other_sp), position=position_dodge(width=0.5), width=0, size=1)+
      geom_text(data=all_p_c, aes(x=Treatment, y=exposure+1, group=other_sp, label=sig),
                position=position_dodge(width=0.4), size=2.5)+
      facet_wrap(~choice, ncol=6)+
      
      scale_y_log10(labels=scales::comma, limits=c(1,1e6))+
      theme_classic(base_size = 18)+
      labs(y=bquote("Daily potential exposures per "~m^2), x="",
           col="Exposing Species")+
      scale_color_manual(values=c("red","black","goldenrod2","brown","hotpink2","gray40"))+
      scale_fill_manual(values=c("red","black","goldenrod2","brown","hotpink2","gray40"))} %>% 
    modify_facet_appearance(strip.background.x.fill = c("red","black","goldenrod2","brown","hotpink2","gray40"),
                            strip.text.x.col = c(rep("white",6)))
  
  
  plot(exp_plot)
  return(exp_plot)
}

# 6b. Facet appearance function needed to change background colors
modify_facet_appearance = function(plot = NULL,
                                     strip.background.x.fill = NULL,
                                     strip.background.y.fill = NULL,
                                     strip.background.x.col = NULL,
                                     strip.background.y.col = NULL,
                                     strip.text.x.col = NULL,
                                     strip.text.y.col = NULL){

   if(is.null(plot)){stop("A ggplot (gg class) needs to be provided!")}

   # Generate gtable object to modify the facet strips:
   g <- ggplot_gtable(ggplot_build(plot))

   # Get the locations of the right and top facets in g:
   stripy <- which(grepl('strip-r|strip-l', g$layout$name)) # account for when strip positions are switched r-l and/or t-b in facet_grid(switch = )
   stripx <- which(grepl('strip-t|strip-b', g$layout$name))

   # Check that the provided value arrays have the same length as strips the plot has:
   lx <- c(length(strip.background.x.fill), length(strip.background.x.col), length(strip.text.x.col))
   if(!all(lx==length(stripx) | lx==0)){stop("The provided vectors with values need to have the same length and the number of facets in the plot!")}
   ly <- c(length(strip.background.y.fill), length(strip.background.y.col), length(strip.text.y.col))
   if(!all(ly==length(stripy) | ly==0)){stop("The provided vectors with values need to have the same length and the number of facets in the plot!")}

   # Change the strips on the y axis:
   for (i in seq_along(stripy)){ # if no strips in the right, the loop will not be executed as seq_along(stripy) will be integer(0)

     # Change strip fill and (border) colour :
     j1 <- which(grepl('strip.background.y', g$grobs[[stripy[i]]]$grobs[[1]]$childrenOrder))
     if(!is.null(strip.background.y.fill[i])){g$grobs[[stripy[i]]]$grobs[[1]]$children[[j1]]$gp$fill <- strip.background.y.fill[i]} # fill
     if(!is.null(strip.background.y.col[i])){g$grobs[[stripy[i]]]$grobs[[1]]$children[[j1]]$gp$col <- strip.background.y.col[i]} # border colour

     # Change color of text:
     j2 <- which(grepl('strip.text.y', g$grobs[[stripy[i]]]$grobs[[1]]$childrenOrder))
     if(!is.null(strip.text.y.col[i])){g$grobs[[stripy[i]]]$grobs[[1]]$children[[j2]]$children[[1]]$gp$col <- strip.text.y.col[i]}

   }

   # Same but for the x axis:
   for (i in seq_along(stripx)){

     # Change strip fill and (border) colour :
     j1 <- which(grepl('strip.background.x', g$grobs[[stripx[i]]]$grobs[[1]]$childrenOrder))
     if(!is.null(strip.background.x.fill[i])){g$grobs[[stripx[i]]]$grobs[[1]]$children[[j1]]$gp$fill <- strip.background.x.fill[i]} # fill
     if(!is.null(strip.background.x.col[i])){g$grobs[[stripx[i]]]$grobs[[1]]$children[[j1]]$gp$col <- strip.background.x.col[i]} # border colour

     # Change color of text:
     j2 <- which(grepl('strip.text.x', g$grobs[[stripx[i]]]$grobs[[1]]$childrenOrder))
     if(!is.null(strip.text.x.col[i])){g$grobs[[stripx[i]]]$grobs[[1]]$children[[j2]]$children[[1]]$gp$col <- strip.text.x.col[i]}

   }

   return(g)
# Note that it returns a gtable object. This can be ploted with plot() or grid::draw.grid().
# patchwork can handle the addition of such gtable to a layout with other ggplot objects,
# but be sure to use patchwork::wrap_ggplot_grob(g) for proper alignment of plots!
# See: https://patchwork.data-imaginist.com/reference/wrap_ggplot_grob.html

}


# 7. significant figures function
sig_fig = function(value, n_figs){
  sf = round(value, n_figs-(1+as.integer(log10(abs(value)))))
  return(sf)
}


# 8. Graphs formatted as network to visualize exposures
make_network_graphs = function(all_e_c_df, metric){
  
   # matrix data
  mepd2 = all_e_c_df %>% 
    dplyr::select(other_sp, Treatment, response, choice) %>% 
    pivot_wider(names_from=other_sp, values_from=response, values_fill = 0) %>%
    mutate(tot = BUFFALO+CATTLE+GIRAFFE+IMPALA+ELEPHANT+ZEBRA) %>%
    mutate_at(vars(BUFFALO:ELEPHANT), funs(./tot)) %>% 
    dplyr::select(Treatment, choice, BUFFALO, CATTLE, GIRAFFE, IMPALA, ELEPHANT,ZEBRA,tot) %>% 
    pivot_longer(c(BUFFALO, CATTLE, GIRAFFE, IMPALA,ELEPHANT,ZEBRA), names_to="exposer", values_to="exposure_prop")
  
  matrix_net_dat = mepd2 %>% 
    filter(Treatment=="Matrix") %>% 
    relocate(exposer, choice, exposure_prop) %>% 
    group_by(exposer, choice) %>%
    summarize_at(vars(exposure_prop), funs(mean(., na.rm=T)))
  
  # water data
  water_net_dat = mepd2 %>% 
    filter(Treatment=="Water") %>% 
    relocate(exposer, choice, exposure_prop) %>% 
    group_by(exposer, choice) %>% 
    summarize_at(vars(exposure_prop), funs(mean(., na.rm=T)))
  
  # create graph
  matrix_net = graph_from_data_frame(matrix_net_dat)
  E(matrix_net)$weight = (matrix_net_dat$exposure_prop)
  E(matrix_net)$color = ifelse(matrix_net_dat$exposer=="CATTLE","red",
                               ifelse(matrix_net_dat$exposer=="BUFFALO","gray20",
                                      ifelse(matrix_net_dat$exposer=="IMPALA","brown",
                                             ifelse(matrix_net_dat$exposer=="GIRAFFE", "goldenrod2",
                                                    ifelse(matrix_net_dat$exposer=="ELEPHANT","gray70","hotpink2")))))
  V(matrix_net)$color=c("gray20","red","gray70","goldenrod2","brown","hotpink2")
  
  water_net = graph_from_data_frame(water_net_dat)
  E(water_net)$weight = (water_net_dat$exposure_prop)
  E(water_net)$color = ifelse(water_net_dat$exposer=="CATTLE","red",
                              ifelse(water_net_dat$exposer=="BUFFALO","gray20",
                                     ifelse(water_net_dat$exposer=="IMPALA","brown",
                                            ifelse(water_net_dat$exposer=="GIRAFFE", "goldenrod2",
                                                   ifelse(water_net_dat$exposer=="ELEPHANT","gray70","hotpink2")))))
  V(water_net)$color=c("gray20","red","gray70","goldenrod2","brown","hotpink2")
  
  
  # remove connections <10%
  matrix_net_pruned = delete.edges(matrix_net, which(E(matrix_net)$weight<0.05))
  water_net_pruned = delete.edges(water_net, which(E(water_net)$weight<0.05))
  
    return(list(water=water_net_pruned, matrix=matrix_net_pruned))
}



# 9. Fit models exploring inter and intra-specific risk
fit_models = function(exp_plot_data, plot_res = F){
  
  # calc risk
  calc_risk = exp_plot_data %>% 
    filter(is.na(exposure)==F) %>% 
    mutate(intra_inter = ifelse(choice==other_sp, "intra", "inter")) %>% 
    group_by(Period, Site, Treatment, choice, intra_inter) %>% 
    summarize_at(vars(exposure), funs(sum(.)))
  
  calc_risk = calc_risk %>% 
    pivot_wider(names_from=intra_inter, values_from=exposure, values_fill = 0) %>% 
    mutate(total = intra+inter) %>% 
    pivot_longer(c(intra, total), names_to="intra_total", values_to="exposure")
  
  # convergence issues if doing intra vs inter rather than intra vs total
  calc_risk$Period = as.factor(calc_risk$Period)
  mod = glmmTMB(exposure ~ Treatment*choice*intra_total + (1|Period)+(1|Site),
                data=calc_risk, family="tweedie")
  
  if(plot_res==T){
    DHARMa::testResiduals(simulateResiduals(mod))
  }
  
  ems = as.data.frame(confint(pairs(emmeans(mod, ~Treatment | intra_total+choice, type="response"),rev=T)))
  pps = as.data.frame(confint(pairs(pairs(emmeans(mod, ~intra_total|Treatment+choice,
                                                  type="response")),
                                    by="choice", rev=F)))
  
  return(list(posthoc_basic = ems, posthoc_pairs = pps, mod=mod))
  
}

# 9b. calc risk function used in #8
calc_risk = function(joined, treatment, activity){
  risked = joined %>% 
    mutate(risk = Ind_Secs*parasite_density,
           new_risk = Ind_Secs*(new_parasites)) %>% 
    filter(Treatment %in% treatment, Activity %in% activity)
  risked = risked %>% 
    group_by(Date, Treatment, Site, choice, new_parasites, parasite_density) %>% 
    summarize_at(vars(Ind_Secs), funs(sum))
  
  risked = risked %>% 
    mutate(new_risk = Ind_Secs*new_parasites,
           risk = Ind_Secs*parasite_density)
  
  return(risked)
}



# 10. Organize the GLMM and marginal means outputs and add a column for the % exposures within water
format_and_plot_intra_total_results = function(output_met, output_lit, output_phy, pct_water){
  all_pps = rbind(output_met$posthoc_pairs,
                  output_lit$posthoc_pairs,
                  output_phy$posthoc_pairs)
 
  all_pps$sharing_metric = c(rep("Metabarcoding", 6), rep("Phylogeny",6), rep("Literature",6))
  
  all_pps = all_pps %>% 
    mutate(sharing_metric = factor(sharing_metric, levels=c("Metabarcoding","Literature","Phylogeny"))) %>% 
    mutate(choice=str_to_title(choice)) %>% 
    mutate(choice=factor(choice, levels=c("Giraffe",
                                          "Buffalo",
                                          "Impala",
                                          "Cattle",
                                          "Elephant",
                                          "Zebra")))
  
  all_pps %>% 
    ggplot(aes(x=choice, y=log(ratio)))+
    geom_point(aes(col=sharing_metric), position=position_dodge(width=0.5),
               size=3)+
    geom_errorbar(aes(col=sharing_metric, ymin=log(lower.CL), ymax=log(upper.CL)),
                  position=position_dodge(width=0.5), width=0, size=1.2)+
    #facet_wrap(~choice)+
    geom_hline(yintercept = 0)+
    theme_classic(base_size=16)+
    theme(legend.position = c(0.8,0.8))+
    labs(x="Focal host species", y="Log-ratio of hotspot effect ln(Inter/Intra)", col="Sharing Metric")+
    scale_color_manual(values=c("black","gray50","gray70"))
  
  
  all_output = rbind(as.data.frame(met_output$posthoc_basic),
                     as.data.frame(lit_output$posthoc_basic),
                     as.data.frame(phy_output$posthoc_basic))
  all_output$metric = c(rep("Metabarcording", 12),
                        rep("Literature", 12),
                        rep("Phylogeny",12))
  
  
  mo = output_met$posthoc_basic %>% 
    as.data.frame() %>% 
    mutate(metric="Metabarcoding")
  lo = output_lit$posthoc_basic %>% 
    as.data.frame() %>% 
    mutate(metric="Literature")
  po = output_phy$posthoc_basic %>% 
    as.data.frame() %>% 
    mutate(metric="Phylogeny")
  ao = rbind(mo, lo, po)
  
  ao = ao %>% 
    mutate_at(vars(metric), funs(ifelse(intra_total=="intra", "intra", .))) %>% 
    group_by(contrast, intra_total, choice, metric) %>% 
    summarize_at(vars(ratio:upper.CL), funs(mean)) %>% 
    mutate(metric2 = paste(intra_total, metric, sep=" - ")) %>% 
    mutate_at(vars(metric2), funs(ifelse(.=="intra - intra", "Intra",.))) %>% 
    mutate(choice=str_to_title(choice), metric2=str_to_title(metric2)) %>% 
    mutate(choice=factor(choice, levels=c("Buffalo","Cattle","Giraffe","Impala",
                                          "Elephant","Zebra"))) %>% 
    mutate(metric2=factor(metric2, levels=c("Intra", "Total - Metabarcoding",
                                            "Total - Literature", "Total - Phylogeny")))
  
  ao %>% 
    ggplot(aes(x=choice, y=ratio))+
    geom_point(aes(col=metric2), position=position_dodge(width=0.5),
               size=2)+
    geom_errorbar(aes(col=metric2, ymin=ratio-SE, ymax=ratio+SE), position=position_dodge(width=0.5),
                  width=0, size=1)+
    theme_classic(base_size=16)+
    labs(x="Focal Species", y="Hotspot Effect (Water/Matrix)", col="Parasite Sharing")+
    scale_color_manual(values=c("red","black","gray50","gray70"))+
    geom_hline(yintercept=1)
  
  
  # back - calculations
  ao$pct_water = calc_ratio_land(ao$ratio,pct_water)
  ao$pct_water_l = calc_ratio_land(ao$ratio-ao$SE, pct_water)
  ao$pct_water_u = calc_ratio_land(ao$ratio+ao$SE, pct_water)
  
  ao %>% 
    ggplot(aes(x=choice, y=pct_water))+
    geom_point(aes(col=metric2), position=position_dodge(width=0.5),
               size=2)+
    geom_errorbar(aes(ymin=pct_water_l, ymax=pct_water_u, col=metric2),
                  position=position_dodge(width=0.5), width=0)+
    theme_classic(base_size=16)+
    scale_color_manual(values=c("red","black","gray50","gray70"))+
    geom_hline(yintercept=pct_water)
  
  return(ao)
}

# 10b. Convert percentage of land
calc_ratio_land = function(ratio, pct){
  r2 = pct*ratio/(100-pct)
  pctw = r2/(r2+1)*100
  return(pctw)
}



# 11. Function that effectively recreates prior analyses (reading everything in from scratch)
# and then reducing cattle parasites to a certain percentage
explore_cattle_psite_reduction = function(pct_psites, pct_water, species){
  
  # modify fecs
  fecs = fecs %>% 
    mutate_at(vars(FEC), funs(ifelse(Species==species, .*pct_psites, .)))
  
  opegg_act = read.csv("data/parasite_risk_at_water.csv")
  
  
  # add FEC measurement
  opegg_act = opegg_act %>% 
    left_join(fecs, by=c("choice"="Species")) %>% 
    mutate(parasite_density = dung_density * FEC) %>% 
    filter(Method == "mpala_median") %>% 
    dplyr::select(-c(dung_density,FEC,Method))
  
  # calculate a 'new parasites' column that adds parasites from other animals
  # weighted by sharing
  op_template = opegg_act %>% 
    filter(Activity=="dgraze") %>% 
    pivot_wider(names_from=choice, values_from = parasite_density, -Ind_Secs) %>% 
    dplyr::select(Period, Date, Treatment, Site, BUFFALO, CATTLE, ELEPHANT, GIRAFFE, IMPALA, ZEBRA)
  
  # repeat for each species
  op_template2 = rbind(op_template, op_template, op_template, op_template, op_template, op_template)
  op_template2$focal_sp = rep(c("Buffalo","Cattle", "Elephant", "Giraffe", "Impala","Zebra"),
                              each=length(op_template$Period))
  
  share_phy = clean_ls(phy_l)
  share_met = clean_ls(met_l)
  share_lit = clean_ls(lit_l)
  
  met_add = calculate_shared_parasites(share_met, op_template2)
  lit_add = calculate_shared_parasites(share_lit, op_template2)
  phy_add = calculate_shared_parasites(share_phy, op_template2)
  
  
  # join tables
  met_join = clean_output(met_add, opegg_act)
  lit_join = clean_output(lit_add, opegg_act)
  phy_join = clean_output(phy_add, opegg_act)
  
  # format
  met_exp_plot_dat = format_join_to_exposure(met_join)
  lit_exp_plot_dat = format_join_to_exposure(lit_join)
  phy_exp_plot_dat = format_join_to_exposure(phy_join)
  
  # fit models
  met_output = fit_models(met_exp_plot_dat)
  phy_output = fit_models(phy_exp_plot_dat)
  lit_output = fit_models(lit_exp_plot_dat)
  
  # means
  met_ems = as.data.frame(emmeans(met_output$mod, ~Treatment | choice+intra_total,
                                  type="response"))
  met_ems$metric = "Metabarcoding"
  phy_ems = as.data.frame(emmeans(phy_output$mod, ~Treatment | choice+intra_total,
                                  type="response"))
  phy_ems$metric = "Phylogeny"
  lit_ems = as.data.frame(emmeans(lit_output$mod, ~Treatment | choice+intra_total,
                                  type="response"))
  lit_ems$metric = "Literature"
  all_ems = rbind(met_ems, phy_ems, lit_ems)
  
  all_pps = rbind(met_output$posthoc_pairs,
                  phy_output$posthoc_pairs,
                  lit_output$posthoc_pairs)
  
  # these are the same?
  comparisons_total_intra = format_and_plot_intra_total_results(met_output, phy_output, lit_output, pct_water)
  comparisons_total_intra$pct_psites = pct_psites
  all_ems$pct_psites = pct_psites
  return(list(comparisons_total_intra = comparisons_total_intra,
              means = all_ems))
  
}


# 12. Function that decreases wildlife parasite density and activity
explore_wildlife_reduction = function(pct_animals, pct_water, species_vec){
  
  # modify fecs (same as modifying dung density here)
  fecs = fecs %>% 
    mutate_at(vars(FEC), funs(ifelse(Species %in% species_vec, .*pct_animals, .)))
  
  opegg_act = read.csv("data/parasite_risk_at_water.csv")
  
  
  # add FEC measurement
  opegg_act = opegg_act %>% 
    left_join(fecs, by=c("choice"="Species")) %>% 
    mutate(parasite_density = dung_density * FEC) %>% 
    filter(Method == "mpala_median") %>% 
    dplyr::select(-c(dung_density,FEC,Method))
  
  # calculate a 'new parasites' column that adds parasites from other animals
  # weighted by sharing
  # also reduce activity by population here
  op_template = opegg_act %>% 
    filter(Activity=="dgraze") %>% 
    mutate_at(vars(Ind_Secs), funs(ifelse(choice %in% species_vec,
                                          .*pct_animals, .))) %>% 
    pivot_wider(names_from=choice, values_from = parasite_density, -Ind_Secs) %>% 
    dplyr::select(Period, Date, Treatment, Site, BUFFALO, CATTLE, ELEPHANT, GIRAFFE, IMPALA, ZEBRA)
  
  # repeat for each species
  op_template2 = rbind(op_template, op_template, op_template, op_template, op_template, op_template)
  op_template2$focal_sp = rep(c("Buffalo","Cattle", "Elephant", "Giraffe", "Impala","Zebra"),
                              each=length(op_template$Period))
  
  share_phy = clean_ls(phy_l)
  share_met = clean_ls(met_l)
  share_lit = clean_ls(lit_l)
  
  met_add = calculate_shared_parasites(share_met, op_template2)
  lit_add = calculate_shared_parasites(share_lit, op_template2)
  phy_add = calculate_shared_parasites(share_phy, op_template2)
  
  
  # join tables
  met_join = clean_output(met_add, opegg_act)
  lit_join = clean_output(lit_add, opegg_act)
  phy_join = clean_output(phy_add, opegg_act)
  
  # format
  met_exp_plot_dat = format_join_to_exposure(met_join)
  lit_exp_plot_dat = format_join_to_exposure(lit_join)
  phy_exp_plot_dat = format_join_to_exposure(phy_join)
  
  # fit models
  met_output = fit_models(met_exp_plot_dat)
  phy_output = fit_models(phy_exp_plot_dat)
  lit_output = fit_models(lit_exp_plot_dat)
  
  # means
  met_ems = as.data.frame(emmeans(met_output$mod, ~Treatment | choice+intra_total,
                                  type="response"))
  met_ems$metric = "Metabarcoding"
  phy_ems = as.data.frame(emmeans(phy_output$mod, ~Treatment | choice+intra_total,
                                  type="response"))
  phy_ems$metric = "Phylogeny"
  lit_ems = as.data.frame(emmeans(lit_output$mod, ~Treatment | choice+intra_total,
                                  type="response"))
  lit_ems$metric = "Literature"
  all_ems = rbind(met_ems, phy_ems, lit_ems)
  
  all_pps = rbind(met_output$posthoc_pairs,
                  phy_output$posthoc_pairs,
                  lit_output$posthoc_pairs)
  
  comparisons_total_intra = format_and_plot_intra_total_results(met_output, phy_output, lit_output, pct_water)
  comparisons_total_intra$pct_animals = pct_animals
  all_ems$pct_psites = pct_animals
  return(list(comparisons_total_intra = comparisons_total_intra,
              means = all_ems))
  
}





