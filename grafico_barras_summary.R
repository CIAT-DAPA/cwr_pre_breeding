### Cargar tablas 

bean_heat        <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress_bean_heat.csv")
carrot_drought   <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress_carrot_drought.csv") 
carrot_heat      <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress_carrot_heat.csv")
eggplant_drought <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress_eggplant_drought.csv")
finger_drought   <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress_finger_millet_drought.csv") 
lentil_drought   <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress_lentil_drought.csv")
potato_drought   <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress_potato_drought.csv")
sunflower_heat   <- read.csv("//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress_sunflower_heat.csv")




db1 <- rbind(bean_heat ,carrot_drought, carrot_heat, eggplant_drought ,finger_drought,lentil_drought, potato_drought, sunflower_heat)
rm (bean_heat ,carrot_drought, carrot_heat, eggplant_drought ,finger_drought,lentil_drought, potato_drought, sunflower_heat) 
saveRDS(db1, "//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress_total.rds")

db1$Yield <- NULL
db1$Percent_yield <- NULL

write.csv(db1,"//dapadfs/Workspace_cluster_9/CWR_pre-breeding/Results/_updated_results/summary_table_crop_stress.csv" )

db1$Crop <- factor(db1$Crop , levels = c("Bean_Heat","Finger_millet_Drought","Potato_Drought","Sunflower_Heat","Lentil_Drought","Eggplant_Drought","Carrot_Drought","Carrot_Heat"))
colours <- c("#1f78b4","#d95f02")
Y <- ggplot(db1, aes(fill=Period, y=(Harvested_area/10^6), x=Crop)) + 
  geom_bar(position="dodge", stat="identity")+
  geom_text(aes(label = Percent_harvested_area, y = Harvested_area/10^6 + 0.15), position = position_dodge(0.9), vjust = 0, color="black", size = 8)+
  scale_fill_manual (values=colours,na.value = "gray")+
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title = element_text(size = 22, face = "bold"),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20)) +
  coord_flip() +
  xlab("Crops and Stress") +
  ylab("Harvested area (million ha)")

ggsave(filename="//dapadfs.cgiarad.org/Workspace_cluster_9/CWR_pre-breeding/Results/harvest.png",plot=Y, width=20, height=10, units='in')





db1$Crop <- factor(db1$Crop , levels = c("Potato_Drought","Eggplant_Drought","Bean_Heat","Finger_millet_Drought","Carrot_Drought","Carrot_Heat","Sunflower_Heat","Lentil_Drought"))
colours <- c("#1f78b4","#d95f02")
Y <- ggplot(db1, aes(fill=Period, y=(Production/10^6), x=Crop)) + 
  geom_bar(position="dodge", stat="identity")+
  # geom_text(aes(label=Percent_yield),vjust=1.6, color="black",
  #           position = position_dodge(0.9), size=8)+
  geom_text(aes(label = Percent_production, y = Production/10^6 + 0.7), position = position_dodge(0.9), vjust = 0, color="black", size = 8)+
  scale_fill_manual (values=colours,na.value = "gray")+
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title = element_text(size = 22, face = "bold"),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20)) +
  coord_flip() +
  xlab("Crops and Stress") +
  ylab("Production (million T)")

ggsave(filename="//dapadfs.cgiarad.org/Workspace_cluster_9/CWR_pre-breeding/Results/production.png",plot=Y, width=20, height=10, units='in')


