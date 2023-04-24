### Analyses for ###

## Cattle aggregations at shared resources
## create potential parasite transmission hotspots for wildlife

# Georgia Titcomb
# georgiatitcomb@gmail.com
# code initiated: 11/21/2021
# last updated: 03/22/2023


##### Setup ####
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(igraph)



# read in data

opegg_act = read.csv("data/parasite_risk_at_water.csv")
all_l = read.csv("data/sharing_matrices.csv")
fecs = read.csv("data/fec_data.csv")

source("scripts/exposure_ratio_functions.R")

# separate sharing metrics
phy_l = all_l %>% filter(type=="phy")
met_l = all_l %>% filter(type=="met")
lit_l = all_l %>% filter(type=="lit")

# add FEC measurement
opegg_act = opegg_act %>% 
  left_join(fecs, by=c("choice"="Species"), multiple="all") %>% 
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

# implement clean function
share_phy = clean_ls(phy_l)
share_met = clean_ls(met_l)
share_lit = clean_ls(lit_l)

# calculate shared parasites
met_add = calculate_shared_parasites(share_met, op_template2)
lit_add = calculate_shared_parasites(share_lit, op_template2)
phy_add = calculate_shared_parasites(share_phy, op_template2)

## Parasite density only (no exposures)
phy_graph_dat = phy_add %>% 
  ungroup() %>% 
  mutate_at(vars(Treatment), funs(recode(Treatment, CONT="Matrix",
                                         WPC="Water")) )%>% 
  pivot_longer(BUFFALO:ZEBRA, names_to="other_sp", values_to="psite") %>% 
  group_by(Treatment, focal_sp, other_sp, Site) %>% 
  summarize_at(vars(psite), funs(mean)) %>% 
  ungroup()

met_graph_dat = met_add %>% 
  ungroup() %>% 
  mutate_at(vars(Treatment), funs(recode(Treatment, CONT="Matrix",
                                         WPC="Water")) )%>% 
  pivot_longer(BUFFALO:ZEBRA, names_to="other_sp", values_to="psite") %>% 
  group_by(Treatment, focal_sp, other_sp, Site) %>% 
  summarize_at(vars(psite), funs(mean)) %>% 
  ungroup()

lit_graph_dat = lit_add %>% 
  ungroup() %>% 
  mutate_at(vars(Treatment), funs(recode(Treatment, CONT="Matrix",
                                         WPC="Water")) )%>% 
  pivot_longer(BUFFALO:ZEBRA, names_to="other_sp", values_to="psite") %>% 
  group_by(Treatment, focal_sp, other_sp, Site) %>% 
  summarize_at(vars(psite), funs(mean)) %>% 
  ungroup()

# add sharing metric
met_graph_dat$sharing_metric="metabarcoding"
lit_graph_dat$sharing_metric="literature"
phy_graph_dat$sharing_metric="phylogeny"

# combine into one df
all_graph_dat = rbind(met_graph_dat, lit_graph_dat, phy_graph_dat)

# format order for plotting
all_graph_dat$focal_sp = factor(all_graph_dat$focal_sp, levels=c("Impala","Giraffe","Buffalo","Zebra",
                                                    "Cattle","Elephant"))
all_graph_dat$other_sp=as.factor(all_graph_dat$other_sp)
levels(all_graph_dat$other_sp)=c("Buffalo","Cattle","Elephant","Giraffe","Impala","Zebra")

# plot all data
all_graph_dat %>% 
  ggplot(aes(x=Treatment, y=(psite+1)))+
  geom_boxplot(aes(fill=other_sp, col=other_sp, alpha=sharing_metric),
               position=position_dodge(width=0.8))+
  geom_point(aes(col=other_sp, alpha=sharing_metric),
             position=position_dodge(width=0.8))+
  scale_fill_manual(values=c("black","orange","gray40","goldenrod","red","orchid"))+
  scale_color_manual(values=c("black","darkorange","gray20","goldenrod4","darkred","orchid4"))+
  facet_wrap(~focal_sp, scales="free")+
  theme_classic(base_size=16)+
  labs(y="Relevant parasite density contributed by each host species", x="",fill="Infected host species",
       color="Infected host species", alpha="Sharing metric")

# using only metabarcoding
all_graph_dat %>% 
  filter(sharing_metric=="metabarcoding") %>% 
  ggplot(aes(x=Treatment, y=psite))+
  geom_boxplot(aes(fill=other_sp, col=other_sp))+
  scale_fill_manual(values=c("black","orange","gray40","goldenrod","red","orchid"))+
  scale_color_manual(values=c("black","darkorange","gray20","goldenrod4","darkred","orchid4"))+
  facet_wrap(~focal_sp, scales="free")+
  theme_classic(base_size=16)+
  labs(y="Parasite density", x="",fill="Source host species",
       color="Source host species")


# clean names and join parasite density and activity tables
met_join = clean_output(met_add, opegg_act)
lit_join = clean_output(lit_add, opegg_act)
phy_join = clean_output(phy_add, opegg_act)


# format
met_exp_plot_dat = format_join_to_exposure(met_join)
lit_exp_plot_dat = format_join_to_exposure(lit_join)
phy_exp_plot_dat = format_join_to_exposure(phy_join)


# Fit models on exposures
# function used:
View(fit_species_specific_exposure_model)

# metabarcoding
b_p_cm = fit_species_specific_exposure_model("BUFFALO", met_exp_plot_dat)
c_p_cm = fit_species_specific_exposure_model("CATTLE", met_exp_plot_dat)
g_p_cm = fit_species_specific_exposure_model("IMPALA", met_exp_plot_dat)
i_p_cm = fit_species_specific_exposure_model("GIRAFFE", met_exp_plot_dat)
e_p_cm = fit_species_specific_exposure_model("ELEPHANT", met_exp_plot_dat)
z_p_cm = fit_species_specific_exposure_model("ZEBRA", met_exp_plot_dat)


# literature
b_p_cl = fit_species_specific_exposure_model("BUFFALO", lit_exp_plot_dat)
c_p_cl = fit_species_specific_exposure_model("CATTLE", lit_exp_plot_dat)
g_p_cl = fit_species_specific_exposure_model("IMPALA", lit_exp_plot_dat)
i_p_cl = fit_species_specific_exposure_model("GIRAFFE", lit_exp_plot_dat)
e_p_cl = fit_species_specific_exposure_model("ELEPHANT", lit_exp_plot_dat)
z_p_cl = fit_species_specific_exposure_model("ZEBRA", lit_exp_plot_dat)

# phylogeny - residuals are imperfect due to small non-zero sharing
b_p_cp = fit_species_specific_exposure_model("BUFFALO", phy_exp_plot_dat)
c_p_cp = fit_species_specific_exposure_model("CATTLE", phy_exp_plot_dat)
g_p_cp = fit_species_specific_exposure_model("IMPALA", phy_exp_plot_dat)
i_p_cp = fit_species_specific_exposure_model("GIRAFFE", phy_exp_plot_dat)
e_p_cp = fit_species_specific_exposure_model("ELEPHANT", phy_exp_plot_dat)
z_p_cp = fit_species_specific_exposure_model("ZEBRA", phy_exp_plot_dat)

all_p_cm = rbind(b_p_cm$pairs, c_p_cm$pairs, g_p_cm$pairs, i_p_cm$pairs)
all_p_cl = rbind(b_p_cl$pairs, c_p_cl$pairs, g_p_cl$pairs, i_p_cl$pairs,
                 e_p_cl$pairs, z_p_cl$pairs)
all_p_cp = rbind(b_p_cp$pairs, c_p_cp$pairs, g_p_cp$pairs, i_p_cp$pairs,
                 e_p_cp$pairs, z_p_cp$pairs)

all_e_cm = rbind(b_p_cm$means, c_p_cm$means, g_p_cm$means, i_p_cm$means,
                 z_p_cm$extra, e_p_cm$extra)
all_e_cl = rbind(b_p_cl$means, c_p_cl$means, g_p_cl$means, i_p_cl$means,
                 z_p_cl$means, e_p_cl$means)
all_e_cp = rbind(b_p_cp$means, c_p_cp$means, g_p_cp$means, i_p_cp$means,
                 z_p_cp$means, e_p_cp$means)

#View(format_and_plot_sp_model_results)
m_plot = format_and_plot_sp_model_results(all_p_cm, all_e_cm, met_exp_plot_dat)
l_plot = format_and_plot_sp_model_results(all_p_cl, all_e_cl, lit_exp_plot_dat)
p_plot = format_and_plot_sp_model_results(all_p_cp, all_e_cp, phy_exp_plot_dat)


# pdf("net_graphs/met_by_species.pdf", width=15, height=5)
plot(m_plot)
# dev.off()


### OUTPUT TABLE 1
options(scipen=0)
met_tab = all_p_cm %>% 
  mutate_at(vars(ratio, t.ratio), funs(as.character(sig_fig(.,4)))) %>% 
  mutate_at(vars(SE, p.value, p.adjusted), funs(as.character(sig_fig(.,3)))) %>% 
  mutate_at(vars(contrast), funs(str_to_title(.))) %>% 
  dplyr::select(contrast, Treatment, ratio, SE, t.ratio,p.adjusted)
lit_tab = all_p_cl %>% 
  mutate_at(vars(ratio, t.ratio), funs(as.character(sig_fig(.,4)))) %>% 
  mutate_at(vars(SE, p.value, p.adjusted), funs(as.character(sig_fig(.,3)))) %>% 
  mutate_at(vars(contrast), funs(str_to_title(.))) %>% 
  dplyr::select(contrast, Treatment, ratio, SE, t.ratio,p.adjusted)
phy_tab = all_p_cp %>% 
  mutate_at(vars(ratio, t.ratio), funs(as.character(sig_fig(.,4)))) %>% 
  mutate_at(vars(SE, p.value, p.adjusted), funs(as.character(sig_fig(.,3)))) %>% 
  mutate_at(vars(contrast), funs(str_to_title(.))) %>% 
  dplyr::select(contrast, Treatment, ratio, SE, t.ratio,p.adjusted)

# library(xlsx)
# write.xlsx(met_tab, "result_tables/species_exposures.xlsx",
#            sheetName = "Metabarcoding", row.names=F)
# write.xlsx(lit_tab, "result_tables/species_exposures.xlsx",
#            sheetName = "Literature", append=T, row.names=F)
# write.xlsx(phy_tab, "result_tables/species_exposures.xlsx",
#            sheetName = "Phylogeny", append=T, row.names=F)


### Network graphs
# for visualization
pms = make_network_graphs(all_e_cm, metric="met")
pls = make_network_graphs(all_e_cl, metric="lit")
pps = make_network_graphs(all_e_cp, metric="phy")

l = layout_in_circle(pms$matrix, order=c("BUFFALO","CATTLE","GIRAFFE","IMPALA","ELEPHANT",
                                                "ZEBRA"))
#pdf("net_graphs/matrix_net_met.pdf", width=8, height=8)
plot.igraph(pms$matrix,
            edge.curved=0.2,
            edge.arrow.size=1.2,
            edge.width=E(pms$matrix)$weight*10,
            edge.color=E(pms$matrix)$color,
            edge.label=sig_fig(E(pms$matrix)$weight,3),
            edge.label.color="black",
            layout=l,
            vertex.color=V(pms$matrix)$color,
            vertex.size=40,
            vertex.label.font=2,
            main="Matrix",
            margin=0.2,
            vertex.label.color="white")
#dev.off()

#pdf("net_graphs/water_net_met.pdf", width=8, height=8)
plot(pms$water,
     edge.curved=0.2,
     edge.arrow.size=1.2,
     edge.arrow.curve=0,
     edge.color=E(pms$water)$color,
     edge.width=E(pms$water)$weight*10,
     edge.label=sig_fig(E(pms$water)$weight,3),
     edge.label.color="black",
     vertex.color=V(pms$water)$color,
     vertex.size=40,
     vertex.label.color="white",
     vertex.label.font=2,
     margin=0.2,
     main="Water",
     layout=l)
#dev.off()



###########################################################
### total risk

## Breakdown of missing data (i.e. failed camera deployments)
met_exp_plot_dat %>%
  filter(is.na(exposure)) %>%
  group_by(Treatment, choice, Period, Site) %>% dplyr::summarize(n=n()) %>% 
  group_by(Treatment, Period, Site) %>% dplyr::summarize(n=n()) %>% 
  group_by(Treatment, Site) %>% dplyr::summarize(n=n())
# 8 periods total for matrix, 5 periods for water (out of 45 total)


# Fit models to compare inter and intraspecific exposures
# (may take a minute to run)
met_output = fit_models(met_exp_plot_dat)
phy_output = fit_models(phy_exp_plot_dat)
lit_output = fit_models(lit_exp_plot_dat)

# combine results
all_pps = rbind(met_output$posthoc_pairs,
                phy_output$posthoc_pairs,
                lit_output$posthoc_pairs)

# format results and add a column to calculate the percentage of exposures near water (0.5%)
comparisons_total_intra = format_and_plot_intra_total_results(met_output, phy_output, lit_output, 0.5)


# plot exposure ratio
comparisons_total_intra %>% 
  mutate(choice = factor(choice, levels=c("Cattle","Buffalo",
                                          "Giraffe","Impala",
                                          "Zebra","Elephant"))) %>% 
  ggplot(aes(x=choice, y=ratio))+
  geom_hline(yintercept = c(10,20,30,40,50,60,70,80,90,100,200,300,400),alpha=0.1)+
  geom_point(aes(col=metric2), position=position_dodge(width=0.5),
             size=2)+
  geom_errorbar(aes(col=metric2, ymin=lower.CL, ymax=upper.CL), position=position_dodge(width=0.5),
                width=0, size=1)+
  theme_classic(base_size=16)+
  labs(x="Focal Species", y="Hotspot Effect (Water/Matrix)", col="Parasite sharing")+
  scale_color_manual(values=c("black","red1","orange2","red4"))+
  geom_vline(xintercept = 4.5)+
  scale_y_log10()+
  geom_hline(yintercept=1)+
  theme(legend.position = c(0.85,0.15))


# plot percent water
comparisons_total_intra %>% 
  ggplot(aes(x=choice, y=pct_water))+
  geom_point(aes(col=metric2), position=position_dodge(width=0.5),
             size=2)+
  geom_errorbar(aes(ymin=pct_water_l, ymax=pct_water_u, col=metric2),
                position=position_dodge(width=0.5), width=0)+
  theme_classic(base_size=16)+
  scale_color_manual(values=c("black","red1","orange2","red4"))+
  geom_hline(yintercept=0.5)


# export table
tab2 = comparisons_total_intra %>% 
  ungroup() %>% 
  arrange(choice, metric2) %>% 
  dplyr::select(choice,metric2, ratio, SE, lower.CL, upper.CL,
                pct_water) %>% 
  mutate_at(vars(ratio, lower.CL, upper.CL),
            funs(as.character(sig_fig(.,4)))) %>% 
  mutate_at(vars(SE, pct_water), funs(as.character(sig_fig(.,3)))) %>% 
  as.data.frame()

# write.xlsx(tab2, "result_tables/intra_total_comparison.xlsx", row.names=F,
#            sheetName = "Sheet1")




#####################################################


### What is the effect of reducing cattle infections?

# vector storing values to consider
pct_consider = c(0.05, seq(from=0.1, to=2, by=0.1))

# calculate very small value (1%!)
comps = explore_cattle_psite_reduction(0.01, 0.5, "CATTLE")
comp0.01 = comps$comparisons_total_intra
mean0.01 = comps$means


# this will take 10 mins to cycle through values stored in vector
for(i in 1:length(pct_consider)){
  comp_temp = explore_cattle_psite_reduction(pct_consider[i], 0.5, "CATTLE")
  comp0.01 = rbind(comp0.01, comp_temp$comparisons_total_intra)
  mean0.01 = rbind(mean0.01, comp_temp$means)
  
}

# preview results
head(comp0.01)


# create plot
# intra hline
comp_intra_mean = comp0.01 %>% 
  filter(choice %in% c("Cattle","Buffalo", "Impala","Giraffe")) %>% 
  filter(metric == "intra") %>% 
  group_by(choice) %>% 
  summarize_at(vars(pct_water), funs(mean)) %>% 
  ungroup()


# plot results
comp0.01 %>% 
  mutate(metric = factor(metric, levels=c("Metabarcoding",
                                  "Literature",
                                  "Phylogeny"))) %>% 
  mutate(metric = recode(metric,
                         "Metabarcoding"="Total - Metabarcoding",
                         "Literature" = "Total - Literature",
                         "Phylogeny" = "Total - Phylogeny")) %>% 
  filter(choice %in% c("Buffalo", "Impala","Giraffe")) %>% 
  filter(metric != "intra") %>% 
  filter(pct_psites != 0) %>% 
  ggplot(aes(x=pct_psites*100, y=pct_water))+
  geom_vline(xintercept = c(0.25,0.5,0.75), alpha=0.05)+
  geom_line(aes(col=metric, linetype=metric))+
  geom_ribbon(aes(col=NULL, fill=metric, ymin=pct_water_l, ymax=pct_water_u),
              alpha=0.1)+
  facet_wrap(~choice, nrow=3)+
  geom_hline(yintercept=0.5,alpha=0.2)+
  # geom_point(data=comp_intra_mean, aes(x=0, y=pct_water), shape=5,
  #            size=2, col="red",stroke=1)+
  scale_color_manual(values=c("red1","orange2","red4"))+
  scale_fill_manual(values=c("red1","orange2","red4"))+
  theme_classic(base_size=16)+
  labs(x="% Cattle parasites", y="% Exposures within 50m water",
       col="Parasite sharing", linetype="Parasite sharing",
       fill="Parasite sharing")



######################
# Wildlife declines


pct_consider = c(0.05, seq(from=0.1, to=2, by=0.2))
baseline = explore_cattle_psite_reduction(1,0.5,"CATTLE")
wild_comp3.5 = explore_wildlife_reduction(0.226, 0.5, c("BUFFALO",
                                                          "IMPALA",
                                                          "GIRAFFE",
                                                          "ZEBRA",
                                                          "ELEPHANT"))
wild_comp8.1 = explore_wildlife_reduction(0.098, 0.5, c("BUFFALO",
                                                        "IMPALA",
                                                        "GIRAFFE",
                                                        "ZEBRA",
                                                        "ELEPHANT"))


# add annotations
baseline$comparisons_total_intra$scenario = "baseline"
wild_comp3.5$comparisons_total_intra$scenario = "wild3.5"
wild_comp8.1$comparisons_total_intra$scenario = "wild8.1"
baseline$means$scenario = "baseline"
wild_comp3.5$means$scenario = "wild3.5"
wild_comp8.1$means$scenario = "wild8.1"

# combine data
compare_wild = rbind(baseline$comparisons_total_intra,
                     wild_comp3.5$comparisons_total_intra,
                     wild_comp8.1$comparisons_total_intra)
compare_wild2 = rbind(baseline$means,
                      wild_comp3.5$means,
                      wild_comp8.1$means)


# plot results
compare_wild %>% 
  filter(choice %in% c("Buffalo","Giraffe","Impala")) %>% 
  mutate(scenario = recode(scenario, "baseline" = "0.8 (OPC baseline)",
                           "wild3.5" = "3.5 (1980)",
                           "wild8.1" = "8.1 (2013)")) %>% 
  ggplot(aes(x=scenario, y=pct_water))+
  geom_point(aes(col=metric2),position=position_dodge(width=0.75),
             size=2)+
  geom_errorbar(aes(ymin=pct_water_l, ymax=pct_water_u, col=metric2),
                position=position_dodge(width=0.75),
                width=0, size=1)+
  facet_wrap(~choice, nrow=3)+
  theme_classic(base_size = 16)+
  scale_color_manual(values=c("black","red1","orange2","red4"))+
  labs(x="Livestock:Wildlife Ratio", y="% Exposures within 50m water",
       col="")+
  #theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))+
  scale_x_discrete(labels=c("0.8\n(OPC 2017)","3.5\n(1980)","8.1\n(2013)"))








