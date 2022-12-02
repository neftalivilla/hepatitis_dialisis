#Prevalent coinfection and associated factors for Hepatitis B, Hepatitis C, and Human Immunodeficiency Virus in patients submitted to renal replacement therapy: a cross-sectional study of 21 dialysis units in the State of Mexico
#Analisis: Neftali E. Antonio Villa, MD
#11 August 2022
#Disclosure: The following code belongs to the authors and should be used for research purpoused. Any issued and clarifications please contact me at neftalivilla@comunidad.unam.mx


setwd("/Users/nefoantonio/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/PROYECTOS/DRA PALOMO/Seroprevalencia Virus Hepatitis")
#---Library Managment-----
library(tidyverse)
library(stringr)
library(haven)
library(readxl)
library(na.tools)
library(parallel)
library(pbmcapply)
library(ggthemes)
library(ggalt)
library(ggpubr)
library(janitor)
library(epiR)
library(glmnet)

#---Dataset Managment-----
base1<- read_excel("Base Evolución VHB y VHC_Efrén_SPP_11 octubre 2019.xlsx", sheet = "BASE COMPLETA")
HEPB<- read_excel("Base Evolución VHB y VHC_Efrén_SPP_11 octubre 2019.xlsx",sheet = "CASOS_HEP_B")
HEPC<- read_excel("Base Evolución VHB y VHC_Efrén_SPP_11 octubre 2019.xlsx",sheet = "CASOS_HEP_C")
base1<- janitor::clean_names(base1)
base1<-base1%>%rename("Folio"=folio_del_paciente)
base1<-base1%>%left_join(HEPB,by="Folio")
base1<-base1%>%left_join(HEPC,by="Folio")

sedes<- read_excel("Datos_Sedes.xlsx")
sedes<-janitor::clean_names(sedes)

#Infeccion Presente
base1$INFECCION_HEP_B<-NULL
base1$INFECCION_HEP_B[base1$h_bs_ag=="POSITIVO" | base1$lg_g_anti_hbc=="POSITIVO" | base1$carga_viral_vhb=="POSITIVO" | ((base1$h_bs_ag=="POSITIVO" | base1$carga_viral_vhb=="POSITIVO") & base1$anti_hbs=="POSITIVO")]<-1
base1$INFECCION_HEP_B[base1$h_bs_ag=="NEGATIVO"]<-0


base1$INFECCION_HEP_C<-NULL
base1$INFECCION_HEP_C[base1$antivirus_hc=="POSITIVO" | base1$arn_virus_hc=="POSITIVO"]<-1
base1$INFECCION_HEP_C[base1$antivirus_hc=="NEGATIVO" & base1$arn_virus_hc=="NEGATIVO"]<-0

base1$INFECCION_VIH<-NULL
base1$INFECCION_VIH[base1$vih=="POSITIVO"]<-1
base1$INFECCION_VIH[base1$vih!="POSITIVO"]<-0

base1$INFECCION_HEP_VIRAL<-NULL
base1$INFECCION_HEP_VIRAL[base1$INFECCION_HEP_B==1 | base1$INFECCION_HEP_C==1]<-1
base1$INFECCION_HEP_VIRAL[base1$INFECCION_HEP_B!=1 & base1$INFECCION_HEP_C!=1]<-0

base1$INFECCION_VIRAL<-NULL
base1$INFECCION_VIRAL[base1$INFECCION_HEP_C==1 | base1$INFECCION_HEP_B==1 | base1$INFECCION_VIH==1]<-1
base1$INFECCION_VIRAL[base1$INFECCION_HEP_C!=1 & base1$INFECCION_HEP_B!=1 & base1$INFECCION_VIH!=1]<-0



#Tiempo de Seguimiento

base1$SEGUIMIENTO_VHB<-NULL
base1$SEGUIMIENTO_VHB<-as.numeric(as.Date(base1$FECHA_DE_DIAGNÓSTICO_HEP_B,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365
base1$SEGUIMIENTO_VHB_CENS<-as.numeric(as.Date(base1$fecha3,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365
base1$SEGUIMIENTO_VHB<-na.tools::na.replace(base1$SEGUIMIENTO_VHB,base1$SEGUIMIENTO_VHB_CENS)
base1$YEAR_DX_VHB<-lubridate::year(as.Date(base1$FECHA_DE_DIAGNÓSTICO_HEP_B,tryFormats = c("%Y-%m-%d")))

base1$SEGUIMIENTO_VHC<-NULL
base1$SEGUIMIENTO_VHC<-as.numeric(as.Date(base1$Fecha_de_DX_HEP_C,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365
base1$SEGUIMIENTO_VHC_CENS<-as.numeric(as.Date(base1$fecha4,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365
base1$SEGUIMIENTO_VHC<-na.tools::na.replace(base1$SEGUIMIENTO_VHC,base1$SEGUIMIENTO_VHC_CENS)
base1$YEAR_DX_VHC<-lubridate::year(as.Date(base1$Fecha_de_DX_HEP_C,tryFormats = c("%Y-%m-%d")))

base1$SEGUIMIENTO_VIH<-NULL
base1$SEGUIMIENTO_VIH<-as.numeric(as.Date(base1$fecha5,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365
base1$SEGUIMIENTO_VIH_CENS<-as.numeric(as.Date(base1$fecha4,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365
base1$SEGUIMIENTO_VIH<-na.tools::na.replace(base1$SEGUIMIENTO_VHC,base1$SEGUIMIENTO_VIH_CENS)

base1$YEAR_DX_VIH<-NULL
base1$YEAR_DX_VIH<-as.numeric(base1$fecha_vih)

base1$SEGUIMIENTO_ANY<-NULL
base1$SEGUIMIENTO_ANY<-as.numeric(as.Date(base1$FECHA_DE_DIAGNÓSTICO_HEP_B,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365
base1$SEGUIMIENTO_ANY<-na.tools::na.replace(base1$SEGUIMIENTO_ANY,as.numeric(as.Date(base1$Fecha_de_DX_HEP_C,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365)
base1$SEGUIMIENTO_ANY<-na.tools::na.replace(base1$SEGUIMIENTO_ANY,as.numeric(as.Date(base1$fecha5,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365)
base1$SEGUIMIENTO_ANY_CENS<-as.numeric(as.Date(base1$fecha4,tryFormats = c("%Y-%m-%d"))-as.Date(base1$fecha_de_inicio_hd,tryFormats = c("%Y-%m-%d")))/365
base1$SEGUIMIENTO_ANY<-na.tools::na.replace(base1$SEGUIMIENTO_ANY,base1$SEGUIMIENTO_ANY_CENS)

base1$YEAR_DX_ANY<-base1$YEAR_DX_VHB
base1$YEAR_DX_ANY<-na.tools::na.replace(base1$YEAR_DX_ANY,base1$YEAR_DX_VHC)
base1$YEAR_DX_ANY<-na.tools::na.replace(base1$YEAR_DX_ANY,base1$YEAR_DX_VIH)

base1$DECADE_INFECCION_VHB<-NULL
base1$DECADE_INFECCION_VHB[base1$YEAR_DX_VHB>=1990 & base1$YEAR_DX_VHB<2000]<-1
base1$DECADE_INFECCION_VHB[base1$YEAR_DX_VHB>=2000 & base1$YEAR_DX_VHB<2010]<-2
base1$DECADE_INFECCION_VHB[base1$YEAR_DX_VHB>=2010 & base1$YEAR_DX_VHB<=2020]<-3

base1$DECADE_INFECCION_VHC<-NULL
base1$DECADE_INFECCION_VHC[base1$YEAR_DX_VHC>=1990 & base1$YEAR_DX_VHC<2000]<-1
base1$DECADE_INFECCION_VHC[base1$YEAR_DX_VHC>=2000 & base1$YEAR_DX_VHC<2010]<-2
base1$DECADE_INFECCION_VHC[base1$YEAR_DX_VHC>=2010 & base1$YEAR_DX_VHC<=2020]<-3

base1$DECADE_INFECCION_VIH<-NULL
base1$DECADE_INFECCION_VIH[base1$YEAR_DX_VIH>=1990 & base1$YEAR_DX_VIH<2000]<-1
base1$DECADE_INFECCION_VIH[base1$YEAR_DX_VIH>=2000 & base1$YEAR_DX_VIH<2010]<-2
base1$DECADE_INFECCION_VIH[base1$YEAR_DX_VIH>=2010 & base1$YEAR_DX_VIH<=2020]<-3

base1$DECADE_INFECCION_ANY<-NULL
base1$DECADE_INFECCION_ANY[base1$DECADE_INFECCION_VHB==1 | base1$DECADE_INFECCION_VHC==1 | base1$DECADE_INFECCION_VIH==1]<-1
base1$DECADE_INFECCION_ANY[base1$DECADE_INFECCION_VHB==2 | base1$DECADE_INFECCION_VHC==2 | base1$DECADE_INFECCION_VIH==2]<-2
base1$DECADE_INFECCION_ANY[base1$DECADE_INFECCION_VHB==3 | base1$DECADE_INFECCION_VHC==3 | base1$DECADE_INFECCION_VIH==3]<-3
base1$DECADE_INFECCION_ANY<-na.tools::na.replace(base1$DECADE_INFECCION_ANY,0)

#Demodialis
base1$modalidad[base1$modalidad=="dpa"]<-"DPA"
base1$modalidad_REC<-NULL
base1$modalidad_REC[base1$modalidad=="DP"]<-1
base1$modalidad_REC[base1$modalidad=="DPA"]<-1
base1$modalidad_REC[base1$modalidad=="DPCA"]<-1
base1$modalidad_REC[base1$modalidad=="DPI"]<-1
base1$modalidad_REC[base1$modalidad=="HD"]<-2

base1$BMI_CAT<-NULL
base1$BMI_CAT[base1$imc<=18.5]<-1
base1$BMI_CAT[base1$imc>18.5 & base1$imc<25]<-2
base1$BMI_CAT[base1$imc>=25 & base1$imc<30]<-3
base1$BMI_CAT[base1$imc>=30]<-4

#Age Categories

base1$EDAD_CAT<-NULL
base1$EDAD_CAT[base1$edad<40]<-1
base1$EDAD_CAT[base1$edad>=40 & base1$edad<65]<-2
base1$EDAD_CAT[base1$edad>=65]<-3

#Sedes 

base2<-base1%>%filter(!is.na(INFECCION_HEP_B))

CASOS_CENTRO_DIALISIS<-base2 %>% group_by(id_centro_de_dialisis) %>%
  summarise(CASOS_PREV_VIRAL_ANY=sum(INFECCION_VIRAL==1),
            CASOS_PREV_HEP_B=sum(INFECCION_HEP_B==1),
            CASOS_PREV_HEP_C=sum(INFECCION_HEP_C==1),
            CASOS_PREV_VIH=sum(INFECCION_VIH==1),
            CASOS_INC_VIRAL_ANY=sum(INFECCION_VIRAL==1 & SEGUIMIENTO_ANY>=0),
            CASOS_INC_HEP_B=sum(INFECCION_HEP_B==1 & SEGUIMIENTO_VHB>=0),
            CASOS_INC_HEP_C=sum(INFECCION_HEP_C==1 & SEGUIMIENTO_VHC>=0),
            CASOS_INC_VIH=sum(INFECCION_VIH==1 & SEGUIMIENTO_VIH>=0))%>%
  mutate(PROP_CASIS=(CASOS_INC_VIRAL_ANY/CASOS_PREV_VIRAL_ANY)*100)%>%
  rename("ID_CENTRO_DIALISIS"="id_centro_de_dialisis")

sedes.2<-sedes%>%rename("ID_CENTRO_DIALISIS"="id_del_centro_de_dialisis")%>%
  left_join(CASOS_CENTRO_DIALISIS,by = "ID_CENTRO_DIALISIS")%>%
  mutate(HEMODIALISIS_FIMMS = coalesce(el_paciente_se_dializa_extramuros_imss,el_paciente_se_hemodializa_extramuros_imss),
         VACUNA_HEP=coalesce(el_personal_de_salud_ha_recibido_la_vacuna_de_hepatitis,el_personal_de_salud_en_dp_ha_recibido_la_vacuna_de_hepatitis),
         MATERIAL_DES=coalesce(utiliza_material_desechable_en_todos_sus_precedimientos,utiliza_material_desechable_en_todos_sus_procedimientos_de_dp),
         SEROLOGIA=coalesce(los_pacientes_cuentan_con_serologia_antes_de_entrar_a_hd,los_pacientes_cuentan_con_serologia_antes_de_ingresar_a_dp),
         FUMIGA=coalesce(fumiga_el_area_con_soluciones_plaguicidas_mensualmente,fumiga_el_area_con_soluciones_plaguicidas_una_vez_al_mes),
         ASEA_SUP=coalesce(asea_las_superficies_del_moviliario_y_equipo_con_detergente,asea_las_superficies_del_moviliario_y_equipo_con_detergente_e_hc),
         DOSIS_VAC=coalesce(cuantas_dosis_de_la_vacuna,cuantas_dosis_recibio_de_la_vacuna),
         AISLAMIENTO=coalesce(el_personal_practica_tecnica_de_aislamiento_y_preven_en_sero,el_personalemplea_metodicamente_tecnicas_de_aislamiento),
         PRACTICAS_HIGIENE_HD=(as.numeric(sedes$utilizan_aguja_arterial_esteril_para_realizar_la_hemodialisis=="SI") +
                                 as.numeric(sedes$utilizan_material_esteril_para_puncionar_las_fistulas=="SI") +
                                 as.numeric(sedes$utiliza_material_desechable_en_todos_sus_precedimientos=="SI") +
                                 as.numeric(sedes$vigila_el_funcionamiento_de_la_maquina_antes_de_la_sesion=="SI") +
                                 as.numeric(sedes$mantiene_seis_horas_la_desinfeccion=="SI") +
                                 as.numeric(sedes$fumiga_el_area_con_soluciones_plaguicidas_mensualmente=="SI") +
                                 as.numeric(sedes$asea_las_superficies_del_moviliario_y_equipo_con_detergente=="SI") +
                                 as.numeric(sedes$utiliza_algun_otro_metodo_de_desinfeccion_de_la_maquina=="SI") +
                                 as.numeric(sedes$tiene_maquinas_para_pacientes_seropositivos_o_desconocidos=="SI") +
                                 as.numeric(sedes$el_personal_practica_tecnica_de_aislamiento_y_preven_en_sero=="SI")+
                                 as.numeric(sedes$los_pacientes_cuentan_con_serologia_antes_de_entrar_a_hd=="SI")+
                                 as.numeric(sedes$el_personal_de_salud_ha_recibido_la_vacuna_de_hepatitis=="SI")),
         PRACTICAS_HIGIENE_DP=(as.numeric(sedes$realiza_lavado_de_manos_antes_y_despues_del_procedimiento_de_dp=="SI") +
                                 as.numeric(sedes$utiliza_material_desechable_en_todos_sus_procedimientos_de_dp=="SI") +
                                 as.numeric(sedes$fumiga_el_area_con_soluciones_plaguicidas_una_vez_al_mes=="SI") +
                                 as.numeric(sedes$asea_las_superficies_del_moviliario_y_equipo_con_detergente_e_hc=="SI") +
                                 as.numeric(sedes$el_personalemplea_metodicamente_tecnicas_de_aislamiento=="SI") +
                                 as.numeric(sedes$los_pacientes_cuentan_con_serologia_antes_de_ingresar_a_dp=="SI") +
                                 as.numeric(sedes$el_personal_de_salud_en_dp_ha_recibido_la_vacuna_de_hepatitis=="SI")),
         NUM_RAC_HIGIENE=coalesce(PRACTICAS_HIGIENE_HD,PRACTICAS_HIGIENE_DP))

sedes.2$CASOS_PREV_VIRAL_ANY[is.na(sedes.2$CASOS_PREV_VIRAL_ANY)]<-0
sedes.2$CASOS_PREV_HEP_B[is.na(sedes.2$CASOS_PREV_HEP_B)]<-0
sedes.2$CASOS_PREV_HEP_C[is.na(sedes.2$CASOS_PREV_HEP_C)]<-0
sedes.2$CASOS_PREV_VIH[is.na(sedes.2$CASOS_PREV_VIH)]<-0

#---Multiple Imputation Analysis#####

base2<-base2%>%mutate(base2, id.2 = rownames(base2))
base3<-base2%>%dplyr::select(id.2,hemoglobina,calcio,sodio_serico,potasio,fosforo,creatinina_serica,urea,albumina_serica)
base2_imp<-mice::mice(base3, m=5, maxit=5,seed = 123)
base2_imp_2<-complete(base2_imp,1)
base2_imp_2<-base2_imp_2%>%na.omit()
base2<-base2%>%left_join(base2_imp_2,by="id.2")%>%
  dplyr::select(-c(hemoglobina.x,calcio.x,sodio_serico.x,potasio.x,fosforo.x,creatinina_serica.x,urea.x,albumina_serica.x))%>%
  rename("hemoglobina"=hemoglobina.y,
         "calcio"=calcio.y,
         "sodio_serico"=sodio_serico.y,
         "potasio"=potasio.y,
         "fosforo"=fosforo.y,
         "creatinina_serica"=creatinina_serica.y,
         "urea"=urea.y,
         "albumina_serica"=albumina_serica.y)

mice::stripplot(base2_imp, pch = 20, cex = 1.2)#ver como quedan las imputaciones junto a cada una de las variables observadas
mice::densityplot(base2_imp)

#Comparisong of both 
summary(base3)
summary(base2_imp_2)

#Missing Values by Percentages
missing.values <- base3 %>%
  gather(key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, isna) %>%
  summarise(num.isna = n()) %>%
  mutate(pct = num.isna / total * 100)

levels <- (missing.values  %>% filter(isna == T) %>%     
             arrange(desc(pct)))$key

percentage.plot <- missing.values %>%
  ggplot() +
  geom_bar(aes(x = reorder(key, desc(pct)), 
               y = pct, fill=isna), 
           stat = 'identity', alpha=0.8) +
  scale_x_discrete(limits = levels) +
  scale_fill_manual(name = "", 
                    values = c('steelblue', 'tomato3'), 
                    labels = c("Present", "Missing")) +
  coord_flip() +
  labs(title = "Percentage of missing values", 
       x = 'Variable', y = "% of missing values")


ggsave(percentage.plot,
       filename = "Percentage_plot.png", bg = "white",
       width = 40, 
       height = 30,
       units=c("cm"),
       dpi = 350,
       limitsize = FALSE) 


#---Estimated Overall Prevalences Rates-----
#Overall Infection

ncas <- table(base1$INFECCION_VIRAL)[2]; npop <- sum(!is.na(base1$INFECCION_HEP_VIRAL))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
         conf.level = 0.95) * 100,2)

#Hepatitis B
ncas <- table(base1$INFECCION_HEP_B)[2]; npop <- sum(!is.na(base1$INFECCION_HEP_B))
tmp <- as.matrix(cbind(ncas, npop))
count.prev2<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
         conf.level = 0.95) * 100,2)

#Hepatitis C
ncas <- table(base1$INFECCION_HEP_C)[2]; npop <- sum(!is.na(base1$INFECCION_HEP_C))
tmp <- as.matrix(cbind(ncas, npop))
count.prev3<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
         conf.level = 0.95) * 100,2)

#HIC
ncas <- table(base1$INFECCION_VIH)[2]; npop <- sum(!is.na(base1$INFECCION_VIH))
tmp <- as.matrix(cbind(ncas, npop))
count.prev4<-paste0(ncas,"/",tmp)[2]
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
         conf.level = 0.95) * 100,2)

prev.fig1.df<-rbind(prev1,prev2,prev3,prev4)
prev.fig1.df$group<-c("Any \nCoinfection","HBV","HCV","HIV")

#---Estimated Prevalences Rates by Sex-----
#Stratification by Sex
#Overall Infection
#Female
ncas <- table(base1[base1$genero=="FEMENINO",]$INFECCION_VIRAL)[2]; npop <- sum(!is.na(base1[base1$genero=="FEMENINO",]$INFECCION_VIRAL))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis B
ncas <- table(base1[base1$genero=="FEMENINO",]$INFECCION_HEP_B)[2]; npop <- sum(!is.na(base1[base1$genero=="FEMENINO",]$INFECCION_HEP_B))
tmp <- as.matrix(cbind(ncas, npop))
count.prev2<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis C
ncas <- table(base1[base1$genero=="FEMENINO",]$INFECCION_HEP_C)[2]; npop <- sum(!is.na(base1[base1$genero=="FEMENINO",]$INFECCION_HEP_C))
tmp <- as.matrix(cbind(ncas, npop))
count.prev3<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#HIC
ncas <- table(base1[base1$genero=="FEMENINO",]$INFECCION_VIH)[2]; npop <- sum(!is.na(base1[base1$genero=="FEMENINO",]$INFECCION_VIH))
tmp <- as.matrix(cbind(ncas, npop))
count.prev4<-paste0(ncas,"/",tmp)[2]
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

prev.fig2.df.sex_fem<-rbind(prev1,prev2,prev3,prev4)
prev.fig2.df.sex_fem$group<-c("Any \nCoinfection","HBV","HCV","HIV")
prev.fig2.df.sex_fem$class<-c("Women")

#Masculino

ncas <- table(base1[base1$genero=="MASCULINO",]$INFECCION_VIRAL)[2]; npop <- sum(!is.na(base1[base1$genero=="MASCULINO",]$INFECCION_VIRAL))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis B
ncas <- table(base1[base1$genero=="MASCULINO",]$INFECCION_HEP_B)[2]; npop <- sum(!is.na(base1[base1$genero=="MASCULINO",]$INFECCION_HEP_B))
tmp <- as.matrix(cbind(ncas, npop))
count.prev2<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis C
ncas <- table(base1[base1$genero=="MASCULINO",]$INFECCION_HEP_C)[2]; npop <- sum(!is.na(base1[base1$genero=="MASCULINO",]$INFECCION_HEP_C))
tmp <- as.matrix(cbind(ncas, npop))
count.prev3<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#HIC
ncas <- table(base1[base1$genero=="MASCULINO",]$INFECCION_VIH)[2]; npop <- sum(!is.na(base1[base1$genero=="MASCULINO",]$INFECCION_VIH))
tmp <- as.matrix(cbind(ncas, npop))
count.prev4<-paste0(ncas,"/",tmp)[2]
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

prev.fig2.df.sex_masc<-rbind(prev1,prev2,prev3,prev4)
prev.fig2.df.sex_masc$group<-c("Any \nCoinfection","HBV","HCV","HIV")
prev.fig2.df.sex_masc$class<-c("Men")
prev.fig2.df.sex<-rbind(prev.fig2.df.sex_fem,prev.fig2.df.sex_masc)

#---Estimated Prevalences Rates by Age of Diagnoses-----
#Stratification by Age
#Overall Infection
#Decade <2000

ncas <- table(base1[base1$DECADE_INFECCION_ANY==1,]$INFECCION_VIRAL)[1]; npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=1,]$INFECCION_VIRAL))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis B
ncas <- table(base1[base1$DECADE_INFECCION_ANY==1,]$INFECCION_HEP_B)[1];  npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=1,]$INFECCION_HEP_B))
tmp <- as.matrix(cbind(ncas, npop))
count.prev2<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis C
ncas <- table(base1[base1$DECADE_INFECCION_ANY==1,]$INFECCION_HEP_C)[2];  ncas<-0;npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=1,]$INFECCION_HEP_C))
tmp <- as.matrix(cbind(ncas, npop))
count.prev3<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#HIC
ncas <- table(base1[base1$DECADE_INFECCION_ANY==1,]$INFECCION_VIH)[2]; ncas<-0; npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=1,]$INFECCION_VIH))
tmp <- as.matrix(cbind(ncas, npop))
count.prev4<-paste0(ncas,"/",tmp)[2]
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

prev.fig4.df.2000<-rbind(prev1,prev2,prev3,prev4)
prev.fig4.df.2000$group<-c("Any \nCoinfection","HBV","HCV","HIV")
prev.fig4.df.2000$class<-c("<2000 years")

#Decade 2000-2010
#Overall Infection

ncas <- table(base1[base1$DECADE_INFECCION_ANY==2,]$INFECCION_VIRAL)[1]; npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=2,]$INFECCION_VIRAL))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis B
ncas <- table(base1[base1$DECADE_INFECCION_ANY==2,]$INFECCION_HEP_B)[2];  npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=2,]$INFECCION_HEP_B))
tmp <- as.matrix(cbind(ncas, npop))
count.prev2<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis C
ncas <- table(base1[base1$DECADE_INFECCION_ANY==2,]$INFECCION_HEP_C)[2];  npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=2,]$INFECCION_HEP_C))
tmp <- as.matrix(cbind(ncas, npop))
count.prev3<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#HIC
ncas <- table(base1[base1$DECADE_INFECCION_ANY==2,]$INFECCION_VIH)[1]; ncas<-0; npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=2,]$INFECCION_VIH))
tmp <- as.matrix(cbind(ncas, npop))
count.prev4<-paste0(ncas,"/",tmp)[2]
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

prev.fig4.df.2000_2<-rbind(prev1,prev2,prev3,prev4)
prev.fig4.df.2000_2$group<-c("Any \nCoinfection","HBV","HCV","HIV")
prev.fig4.df.2000_2$class<-c("2000 to 2010 decade")


#≥2010 years
#Overall Infection

ncas <- table(base1[base1$DECADE_INFECCION_ANY==3,]$INFECCION_VIRAL)[1]; npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=3,]$INFECCION_VIRAL))
tmp <- as.matrix(cbind(ncas, npop))
count.prev1<-paste0(ncas,"/",tmp)[2]
prev1<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis B
ncas <- table(base1[base1$DECADE_INFECCION_ANY==3,]$INFECCION_HEP_B)[2];  npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=3,]$INFECCION_HEP_B))
tmp <- as.matrix(cbind(ncas, npop))
count.prev2<-paste0(ncas,"/",tmp)[2]
prev2<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#Hepatitis C
ncas <- table(base1[base1$DECADE_INFECCION_ANY==3,]$INFECCION_HEP_C)[2];  npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=3,]$INFECCION_HEP_C))
tmp <- as.matrix(cbind(ncas, npop))
count.prev3<-paste0(ncas,"/",tmp)[2]
prev3<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

#HIC
ncas <- table(base1[base1$DECADE_INFECCION_ANY==3,]$INFECCION_VIH)[2]; npop <- sum(!is.na(base1[base1$DECADE_INFECCION_ANY!=3,]$INFECCION_VIH))
tmp <- as.matrix(cbind(ncas, npop))
count.prev4<-paste0(ncas,"/",tmp)[2]
prev4<-round(epi.conf(tmp, ctype = "prevalence", method = "wilson", N = 10000, design = 1, 
                      conf.level = 0.95) * 100,2)

prev.fig4.df.2000_3<-rbind(prev1,prev2,prev3,prev4)
prev.fig4.df.2000_3$group<-c("Any \nCoinfection","HBV","HCV","HIV")
prev.fig4.df.2000_3$class<-c(">2010 years")
prev.fig4.df.year<-rbind(prev.fig4.df.2000,prev.fig4.df.2000_2,prev.fig4.df.2000_3)
prev.fig4.df.year$class<-factor(prev.fig4.df.year$class,levels = c("<2000 years","2000 to 2010 decade",">2010 years"))

#---Merging Prevalence  (Figure 1)-----

Figure1A<-ggplot(prev.fig1.df, aes(x=group, y=est, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  labs(fill="Viral Infection")+
  xlab("")+
  ylab("Prevalence, (%)")+
  ggsci::scale_fill_jama()+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 1.3),
    position = position_dodge(0.9),
    vjust = 0)+
  scale_y_continuous(limits = c(-.1,9))


Figure1B.1<-ggplot(prev.fig2.df.sex,aes(x=group, y=est, group=class,fill=class))+
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  labs(fill="Sex")+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Dark2")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 1.8),
    position = position_dodge(0.9),
    vjust = 0,size=2.5)+
  scale_y_continuous(limits = c(-.1,9))
  
Figure1B.2<-ggplot(prev.fig4.df.year,aes(x=group, y=est, group=class,fill=class))+
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  theme_pubclean()+
  labs(fill="Decade of Diagnosis")+
  xlab("")+
  ylab("Prevalence, (%)")+
  scale_fill_brewer(palette = "Set1")+ geom_text(
    aes(label = paste0(est,"%","\n","(",lower,"-",upper,")"), y = est + 1.8),
    position = position_dodge(1),
    vjust = 0,size=2.2)+
  scale_y_continuous(limits = c(-.1,9))+
  guides(fill = guide_legend(ncol = 3, nrow = 1, byrow = TRUE))

Figure1A_B<-ggarrange(Figure1B.1,Figure1B.2,ncol = 1,nrow = 2,labels = c("B","C"),common.legend = F)
Figure1<-ggarrange(Figure1A,Figure1A_B,ncol = 2,nrow = 1,common.legend = F)

ggsave(Figure1,
       filename = "Figure1.pdf", 
       width = 40, 
       height = 20,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)


#---Characterization of the center of hemodyalisis (Table 2)-----

#Table Creation
columns <- c('Parameter',"All-Centers (n=21)")
#Modalidad
mod <- table(sedes.2$modalidad2,useNA = "always")[c(1,2)]
mod.prop <- round(prop.table(table(sedes.2$modalidad2,useNA = "always")),4)[c(1,2)]*100
MODALIDAD<-`names<-`(data.frame(c("Peritoneal Dyalisis (%)","Hemodyalisis (%)"),
                                matrix(c(paste(mod,paste0('(',mod.prop,')'))),ncol = 1)),
                     columns)

#Modalidad
hemo.FIMMS <- table(sedes.2$HEMODIALISIS_FIMMS,useNA = "always")[c(2)]
hemo.FIMMS.prop <- round(prop.table(table(sedes.2$HEMODIALISIS_FIMMS,useNA = "always")),4)[c(2)]*100
FIMMS<-`names<-`(data.frame(c("Ambulant RRT (%)"),
                            matrix(c(paste(hemo.FIMMS,paste0('(',hemo.FIMMS.prop,')'))),ncol = 1)),
                 columns)

#Material Desinfeccion
mat <- table(sedes.2$MATERIAL_DES,useNA = "always")[c(2)]
mat.prop <- round(prop.table(table(sedes.2$MATERIAL_DES,useNA = "always")),4)[c(2)]*100
MATERIAL_DESIN<-`names<-`(data.frame(c("Use of Disposable Material (%)"),
                                     matrix(c(paste(mat,paste0('(',mat.prop,')'))),ncol = 1)),
                          columns)

#Serologia
sero <- table(sedes.2$SEROLOGIA,useNA = "always")[c(2)]
sero.prop <- round(prop.table(table(sedes.2$SEROLOGIA,useNA = "always")),4)[c(2)]*100
SEROLOGIA<-`names<-`(data.frame(c("Serology Before RRT(%)"),
                                matrix(c(paste(sero,paste0('(',sero.prop,')'))),ncol = 1)),
                     columns)

#Fumigacion
fumiga <- table(sedes.2$FUMIGA,useNA = "always")[c(2)]
fumiga.prop <- round(prop.table(table(sedes.2$FUMIGA,useNA = "always")),4)[c(2)]*100
FUMIGACION<-`names<-`(data.frame(c("Montly Fumigation (%)"),
                                 matrix(c(paste(fumiga,paste0('(',fumiga.prop,')'))),ncol = 1)),
                      columns)

#Clean Area
asea <- table(sedes.2$ASEA_SUP,useNA = "always")[c(2)]
asea.prop <- round(prop.table(table(sedes.2$ASEA_SUP,useNA = "always")),4)[c(2)]*100
ASEO<-`names<-`(data.frame(c("Periodic Cleaning (%)"),
                           matrix(c(paste(asea,paste0('(',asea.prop,')'))),ncol = 1)),
                columns)

#Vaccination
vac <- table(sedes.2$VACUNA_HEP,useNA = "always")[c(2)]
vac.prop <- round(prop.table(table(sedes.2$VACUNA_HEP,useNA = "always")),4)[c(2)]*100
VACUNA<-`names<-`(data.frame(c("Hepatitis Vacctination (%)"),
                             matrix(c(paste(vac,paste0('(',vac.prop,')'))),ncol = 1)),
                  columns)

#Vaccination Number
vac.num <- table(sedes.2$DOSIS_VAC,useNA = "always")[c(2,3,4,5)]
vac.num.prop <- round(prop.table(table(sedes.2$DOSIS_VAC,useNA = "always")),4)[c(2,3,4,5)]*100
VACUNA_NUM<-`names<-`(data.frame(c("1 (%)", "2 (%)", "3 (%)", "4 (%)"),
                                 matrix(c(paste(vac.num,paste0('(',vac.num.prop,')'))),ncol = 1)),
                      columns)

#Aislamiento
aisla <- table(sedes.2$AISLAMIENTO,useNA = "always")[c(1)]
aisla.prop <- round(prop.table(table(sedes.2$AISLAMIENTO,useNA = "always")),4)[c(1)]*100
AISLAMIENTO<-`names<-`(data.frame(c("Isolation Techniques (%)"),
                                  matrix(c(paste(aisla,paste0('(',aisla.prop,')'))),ncol = 1)),
                       columns)

Table6.1<-rbind(MODALIDAD,MATERIAL_DESIN,SEROLOGIA,FUMIGACION,ASEO,AISLAMIENTO,VACUNA,VACUNA_NUM)
Table6_Flex<-flextable::align(flextable::flextable(Table6.1,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table6_Flex,path="Table_6.docx")


#---Risk Factors Associated with Viral Hepatitis (Figure 3)----

#Modelos de Familiares y Conductas de los Pacientes
base2<-base1%>%filter(!is.na(INFECCION_HEP_B))

familiares<-base2%>%
  group_by(id_centro_de_dialisis)%>%
  summarise(FAMILIARES_HEPATITIS=sum(x15_existe_alguna_persona_que_tenga_hepatitis_en_su_familia=="SI", na.rm = T),
            TIEMPO_FAM_HEPATITIS=mean(x16_tiempo2,na.rm=T),
            RECIBIR_TX=sum(x17_recibio_tratamiento=="SI", na.rm = T),
            OTROS_FAM_VAC_HEPATITIS_B=sum(x18_sus_fam_y_cuidadores_cuentan_con_la_vacuna_de_hepatitis_b=="SI", na.rm = T),
            PADECIO_HEPATITIS=sum(x19_pedecio_hepatitis=="SI", na.rm = T),
            TIEMPO_PADECIO_HEPATITIS=mean(x22_tiempo=="SI", na.rm = T),
            TIPO_HEPATITIS_A=sum(x20_que_tipo_de_hepatitis=="A", na.rm = T),
            TIPO_HEPATITIS_B=sum(x20_que_tipo_de_hepatitis=="B", na.rm = T),
            TIPO_HEPATITIS_C=sum(x20_que_tipo_de_hepatitis=="C", na.rm = T),
            BIOPSIA_HEPATICA=sum(x23_le_realizaron_biopsia_hepatica=="SI", na.rm = T),
            TENIA_VAC_ANTES_DX=sum(x24_tenia_la_vacuna_de_hepatits_antes_del_diagnostico=="SI", na.rm = T),
            TRANSFUCIONES=sum(x26_ha_recibido_transfuciones=="SI", na.rm = T),
            UTILIZA_DROGA=sum(x29_utiliza_algun_tipo_de_droga=="SI", na.rm = T),
            UTILIZA_ESTIMULANTES=sum(x31_utiliza_estimulantes_inyectados=="SI", na.rm = T),
            COMPARTE_AGUJAS=sum(x32_comparte_agujas=="SI", na.rm = T),
            UTILIZA_ACUPUNTURA=sum(x33_utiliza_acupuntura=="SI", na.rm = T),
            TATUAJES=sum(x36_tiene_tatuajes=="SI", na.rm = T),
            PERFORACIONES=sum(x38_perforciones=="SI", na.rm = T),
            CEPILLO_DENTAL=sum(x40_comparte_cepillo_dental=="SI", na.rm = T),
            COMPARTE_AFEITAR=sum(x42_comparte_instrumentos_para_afeitar=="SI", na.rm = T),
            PRUEBAS_VIH=sum(x44_le_han_realizado_pruebas_de_vih=="SI", na.rm = T),
            PAREJAS_SEX=mean(x49_numero_de_parejas_sexuales, na.rm = T))%>%
  rename("ID_CENTRO_DIALISIS"="id_centro_de_dialisis")

sedes.3<-sedes.2%>%left_join(familiares,by="ID_CENTRO_DIALISIS")


#Poisson Regression models 
#Caracteristicas de la Clinica

sedes.3$FUMIGA_ASEA_DESINFECTA<-as.numeric(sedes.3$FUMIGA=="SI")+as.numeric(sedes.3$ASEA_SUP=="SI")+as.numeric(sedes.3$MATERIAL_DES=="SI")
sedes.3$DIALISIS_PERI[sedes.3$modalidad2=="DP"]<-1
sedes.3$DIALISIS_PERI[sedes.3$modalidad2=="HD"]<-0

mod1<-glm(CASOS_PREV_VIRAL_ANY~DIALISIS_PERI+HEMODIALISIS_FIMMS+SEROLOGIA+FUMIGA_ASEA_DESINFECTA+DOSIS_VAC+offset(log(total_de_pacientes)),family = poisson(link = "log"),data = sedes.3)
summary(mod1)
jtools::summ(mod1,exp=T)
car::vif(mod1)

#Figura Familares
Figure3A<-jtools::plot_summs(mod1, coefs = 
                               c("Peritoneal Dialysis Unit"="DIALISIS_PERI",
                                 "Surrogated Hemodialysis Unit"="HEMODIALISIS_FIMMSSI",
                                 "Hepatitis Serologic Testing \nPrior Starting RRT"="SEROLOGIASI",
                                 "Frequently Fumigation, Cleaning \nand Disinfection"="FUMIGA_ASEA_DESINFECTA",
                                 "Number of HBV Vaccines"="DOSIS_VAC"),colors = c("#00b4d8") , exp=T) + 
  xlab("Adjusted Poisson Regression Models: \nIRR, 95%CI")+
  scale_x_log10()+ylab("")+ggtitle("Treatment Facility Associated Factors")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  ylab("")+
  labs(caption = "")+
  theme(legend.position = "top")

#Modelo Familiares
base2$familiar_hepatitis<-NULL
base2$familiar_hepatitis[base2$x15_existe_alguna_persona_que_tenga_hepatitis_en_su_familia=="SI" | base2$x16_tiempo2>0 | base2$x17_recibio_tratamiento=="SI"]<-1
base2$familiar_hepatitis<-na.tools::na.replace(base2$familiar_hepatitis,0)

base2$PERFOS<-NULL
base2$PERFOS[base2$x36_tiene_tatuajes=="SI" | base2$x38_perforciones=="SI" | base2$x33_utiliza_acupuntura=="SI"]<-1
base2$PERFOS<-na.tools::na.replace(base2$PERFOS,0)

mod2<-glm(INFECCION_VIRAL~familiar_hepatitis,
          family = "binomial",data = base2)


summary(mod2)
jtools::summ(mod2,exp=T)
car::vif(mod2)


Figure3<-ggarrange(Figure3A,ncol = 1,nrow = 1,labels = LETTERS[1:2])

ggsave(Figure3,
       filename = "Figure3.pdf", 
       width = 20, height = 15,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)


 #---Characteristics of studies population (Table 1, Supplementary Table 1)-----

base2<-base1%>%filter(!is.na(INFECCION_HEP_B))
columns <- c('Parameter',"All-Population (n=1,304)")

#Sexo

sexo <- table(base2$genero,useNA = "always")[c(2)]
sexoprop <- round(prop.table(table(base2$genero,useNA = "always")),4)[c(2)]*100
Sexo<-`names<-`(data.frame("Men (%)",
                           matrix(c(paste(sexo,paste0('(',sexoprop,')'))),ncol = 1)),
                columns)

#Edad

Edad.media<-round(mean(base1$edad,na.rm = T),2)
Edad.sd<-round(sd(base1$edad,na.rm = T),2)
Edad<-`names<-`(data.frame(matrix(c("Age (years)",paste(Edad.media,"(","±",Edad.sd,")")),ncol = 2,byrow = T)),columns)

#Modalidad TSR

hemodialisis <- table(base2$modalidad_REC,useNA = "always")[c(1:2)]
hemodialisis.prop <- round(prop.table(table(base2$modalidad_REC,useNA = "always")),4)[c(1:2)]*100
Hemodial<-`names<-`(data.frame(c("Peritoneal Dialysis (%)","Hemodialysis (%)"),
                               matrix(c(paste(hemodialisis,paste0('(',hemodialisis.prop,')'))),ncol = 1,byrow = F)),
                    columns)

#Glomerulopatias
glomeru <- table(base2$glomerulopatias,useNA = "always")[c(2)]
glomeru.prop <- round(prop.table(table(base2$glomerulopatias,useNA = "always")),4)[c(2)]*100
Glomerulopatias<-`names<-`(data.frame("Glomerulopathies (%)",
                                      matrix(c(paste(glomeru,paste0('(',glomeru.prop,')'))),ncol = 1)),
                           columns)

#Nefropatia Hipertensiva
nefrop.has <- table(base2$nefropatia_hipertensiva,useNA = "always")[c(2)]
nefrop.has.prop <- round(prop.table(table(base2$nefropatia_hipertensiva,useNA = "always")),4)[c(2)]*100
Nefro.HAS<-`names<-`(data.frame("Hypertensive Nefropathy (%)",
                                matrix(c(paste(nefrop.has,paste0('(',nefrop.has.prop,')'))),ncol = 1)),
                     columns)

#Lupus
lupus <- table(base2$lupus_eritematoso_sistemico,useNA = "always")[c(2)]
lupus.prop <- round(prop.table(table(base2$lupus_eritematoso_sistemico,useNA = "always")),4)[c(2)]*100
Lupus<-`names<-`(data.frame("Lupus (%)",
                            matrix(c(paste(lupus,paste0('(',lupus.prop,')'))),ncol = 1)),
                 columns)

#Lupus
hipoplasia <- table(base2$hipoplasia_renal,useNA = "always")[c(2)]
hipoplasia.prop <- round(prop.table(table(base2$hipoplasia_renal,useNA = "always")),4)[c(2)]*100
Hipoplasia.renal<-`names<-`(data.frame("Hypoplastic Nepropathy (%)",
                                       matrix(c(paste(hipoplasia,paste0('(',hipoplasia.prop,')'))),ncol = 1)),
                            columns)

#Diabetes Nefropatia
diab <- table(base2$diabetes_melitus,useNA = "always")[c(2)]
diab.prop <- round(prop.table(table(base2$diabetes_melitus,useNA = "always")),4)[c(2)]*100
Diabetes.Nefro<-`names<-`(data.frame("Diabetic Nepropathy(%)",
                                     matrix(c(paste(diab,paste0('(',diab.prop,')'))),ncol = 1)),
                          columns)

#Riñon Poliquistico
poliqu <- table(base2$rinon_poliquistico,useNA = "always")[c(2)]
poliqu.prop <- round(prop.table(table(base2$rinon_poliquistico,useNA = "always")),4)[c(2)]*100
Riñon.Poli<-`names<-`(data.frame("Polycystic Kidney Disease (%)",
                                 matrix(c(paste(poliqu,paste0('(',poliqu.prop,')'))),ncol = 1)),
                      columns)

#Uropatia Obstructiva
obstructiva <- table(base2$uropatia_obstructuva_reflujo,useNA = "always")[c(2)]
obstructiva.prop <- round(prop.table(table(base2$uropatia_obstructuva_reflujo,useNA = "always")),4)[c(2)]*100
Uropatia.Obs<-`names<-`(data.frame("Obstructive Uropathy (%)",
                                   matrix(c(paste(obstructiva,paste0('(',obstructiva.prop,')'))),ncol = 1)),
                        columns)

#Nefropatia Uratos
nefro.uro <- table(base2$nefropatia_por_uratos,useNA = "always")[c(2)]
nefro.uro.prop <- round(prop.table(table(base2$nefropatia_por_uratos,useNA = "always")),4)[c(2)]*100
Nefropatia.Uratos<-`names<-`(data.frame("Urate Nephropathy (%)",
                                        matrix(c(paste(nefro.uro,paste0('(',nefro.uro.prop,')'))),ncol = 1)),
                             columns)

#Desconocida
desconocida <- table(base2$desconocido,useNA = "always")[c(2)]
desconocida.prop <- round(prop.table(table(base2$desconocido,useNA = "always")),4)[c(2)]*100
Nefro.desconocida<-`names<-`(data.frame("Unclassified Nephropathy (%)",
                                        matrix(c(paste(desconocida,paste0('(',desconocida.prop,')'))),ncol = 1)),
                             columns)

#Tabaquismo
tabaq <- table(base2$tabaquismo_actual,useNA = "always")[c(2)]
tabaq.prop <- round(prop.table(table(base2$tabaquismo_actual,useNA = "always")),4)[c(2)]*100
Tabaquismo<-`names<-`(data.frame("Current Smoking (%)",
                                 matrix(c(paste(tabaq,paste0('(',tabaq.prop,')'))),ncol = 1)),
                      columns)

#Alcoholismo
alcoholismo <- table(base2$alcoholismo,useNA = "always")[c(2)]
alcoholismo.prop <- round(prop.table(table(base2$alcoholismo,useNA = "always")),4)[c(2)]*100
Alcoholismo<-`names<-`(data.frame("Alcoholism (%)",
                                  matrix(c(paste(alcoholismo,paste0('(',alcoholismo.prop,')'))),ncol = 1)),
                       columns)

#Diabetes
diabetes <- table(base2$dm,useNA = "always")[c(2)]
diabetes.prop <- round(prop.table(table(base2$dm,useNA = "always")),4)[c(2)]*100
Diabetes<-`names<-`(data.frame("Diabetes (%)",
                               matrix(c(paste(diabetes,paste0('(',diabetes.prop,')'))),ncol = 1)),
                    columns)

#Arterial Hypertension
hiper.artr <- table(base2$hta,useNA = "always")[c(2)]
hiper.artr.prop <- round(prop.table(table(base2$hta,useNA = "always")),4)[c(2)]*100
HAS<-`names<-`(data.frame("Arterial Hypertension (%)",
                          matrix(c(paste(hiper.artr,paste0('(',hiper.artr.prop,')'))),ncol = 1)),
               columns)

#Previous Cardiovascular Disease
cardivas <- table(base2$enf_cardiovascular,useNA = "always")[c(2)]
cardivas.prop <- round(prop.table(table(base2$enf_cardiovascular,useNA = "always")),4)[c(2)]*100
CVD<-`names<-`(data.frame("Previous CVD (%)",
                          matrix(c(paste(cardivas,paste0('(',cardivas.prop,')'))),ncol = 1)),
               columns)

#Retinopatia
retino <- table(base2$retinopatia,useNA = "always")[c(2)]
retino.prop <- round(prop.table(table(base2$retinopatia,useNA = "always")),4)[c(2)]*100
Retinopatia<-`names<-`(data.frame("Retinopathy (%)",
                                  matrix(c(paste(retino,paste0('(',retino.prop,')'))),ncol = 1)),
                       columns)

#Alteracion Neurologica
neurologica <- table(base2$alteracion_neurologica,useNA = "always")[c(2)]
neurologica.prop <- round(prop.table(table(base2$alteracion_neurologica,useNA = "always")),4)[c(2)]*100
Neurologic<-`names<-`(data.frame("Neurological Impairment (%)",
                                 matrix(c(paste(neurologica,paste0('(',neurologica.prop,')'))),ncol = 1)),
                      columns)

#Anemia
anemia <- table(base2$anemia,useNA = "always")[c(2)]
anemia.prop <- round(prop.table(table(base2$anemia,useNA = "always")),4)[c(2)]*100
Anemia<-`names<-`(data.frame("Anemia (%)",
                             matrix(c(paste(anemia,paste0('(',anemia.prop,')'))),ncol = 1)),
                  columns)

#Eritropoyetina
eritro <- table(base2$utiliza_eritropoyetina,useNA = "always")[c(2)]
eritro.prop <- round(prop.table(table(base2$utiliza_eritropoyetina,useNA = "always")),4)[c(2)]*100
EPO<-`names<-`(data.frame("Erythropoietin use (%)",
                          matrix(c(paste(eritro,paste0('(',eritro.prop,')'))),ncol = 1)),
               columns)

#Calcioantagonistas
calcio.bloq <- table(base2$calcioantagonista,useNA = "always")[c(2)]
calcio.bloq.prop <- round(prop.table(table(base2$calcioantagonista,useNA = "always")),4)[c(2)]*100
Calcioantago<-`names<-`(data.frame("CCBs use (%)",
                                   matrix(c(paste(calcio.bloq,paste0('(',calcio.bloq.prop,')'))),ncol = 1)),
                        columns)

#Betabloqueadores
beta.bloq <- table(base2$betabloqueador,useNA = "always")[c(2)]
beta.bloq.prop <- round(prop.table(table(base2$betabloqueador,useNA = "always")),4)[c(2)]*100
Betabloq<-`names<-`(data.frame("Beta-blockers use (%)",
                               matrix(c(paste(beta.bloq,paste0('(',beta.bloq.prop,')'))),ncol = 1)),
                    columns)

#Tiazide Use
tiazida.bloq <- table(base2$diuretico_tiazida,useNA = "always")[c(2)]
tiazida.bloq.prop <- round(prop.table(table(base2$diuretico_tiazida,useNA = "always")),4)[c(2)]*100
Tiazidas<-`names<-`(data.frame("Thiazides use (%)",
                               matrix(c(paste(tiazida.bloq,paste0('(',tiazida.bloq.prop,')'))),ncol = 1)),
                    columns)

#Diureticos de Asa
diuretico.asa <- table(base2$diuretico_de_asa,useNA = "always")[c(2)]
diuretico.asa.prop <- round(prop.table(table(base2$diuretico_de_asa,useNA = "always")),4)[c(2)]*100
Diuretico.Asa<-`names<-`(data.frame("Loop Diuretics use (%)",
                                    matrix(c(paste(diuretico.asa,paste0('(',diuretico.asa.prop,')'))),ncol = 1)),
                         columns)

#IECA
eca.bloq <- table(base2$eca,useNA = "always")[c(2)]
eca.bloq.prop <- round(prop.table(table(base2$eca,useNA = "always")),4)[c(2)]*100
IECA<-`names<-`(data.frame("ACEI use (%)",
                           matrix(c(paste(eca.bloq,paste0('(',eca.bloq.prop,')'))),ncol = 1)),
                columns)

#Alfa Bloq
alfa.bloq <- table(base2$alfabloqueador,useNA = "always")[c(2)]
alfa.bloq.prop <- round(prop.table(table(base2$alfabloqueador,useNA = "always")),4)[c(2)]*100
Alfa.Bloqueador<-`names<-`(data.frame("Alpha Blockers use (%)",
                                      matrix(c(paste(alfa.bloq,paste0('(',alfa.bloq.prop,')'))),ncol = 1)),
                           columns)

#ARA-II Bloq
ARAII.bloq <- table(base2$ara_ii,useNA = "always")[c(2)]
ARAII.prop <- round(prop.table(table(base2$ara_ii,useNA = "always")),4)[c(2)]*100
ARAII.bloqueador<-`names<-`(data.frame("ARBs use (%)",
                                       matrix(c(paste(ARAII.bloq,paste0('(',ARAII.prop,')'))),ncol = 1)),
                            columns)

#Hierro Suplementario
hierro <- table(base2$hierro,useNA = "always")[c(2)]
hierro.prop <- round(prop.table(table(base2$hierro,useNA = "always")),4)[c(2)]*100
Hierro.Sup<-`names<-`(data.frame("Iron Supplement Intake (%)",
                                 matrix(c(paste(hierro,paste0('(',hierro.prop,')'))),ncol = 1)),
                      columns)

#Uso Hipoglucemiantes
hipoglu <- table(base2$hipoglucemiantes,useNA = "always")[c(2)]
hipoglu.prop <- round(prop.table(table(base2$hipoglucemiantes,useNA = "always")),4)[c(2)]*100
Hipoglucemiantes<-`names<-`(data.frame("Hypoglycemiant use (%)",
                                       matrix(c(paste(hipoglu,paste0('(',hipoglu.prop,')'))),ncol = 1)),
                            columns)

#Multivitamins
multivit <- table(base2$multivitaminas,useNA = "always")[c(2)]
multivit.prop <- round(prop.table(table(base2$multivitaminas,useNA = "always")),4)[c(2)]*100
Multivitamins<-`names<-`(data.frame("Multivitamins use (%)",
                                    matrix(c(paste(multivit,paste0('(',multivit.prop,')'))),ncol = 1)),
                         columns)

#Peso

num1<-c(paste(round(median(base1$peso_actual,na.rm = T ),2),
              paste0('(',round(quantile(base1$peso_actual,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$peso_actual,na.rm = T,probs = c(0.75)),2),')')))

PESO<-`names<-`(data.frame(matrix(c("Weight (kg)",num1),ncol = 2)),columns)

#Talla

Talla.media<-round(mean(base1$talla,na.rm = T),2)
Talla.sd<-round(sd(base1$talla,na.rm = T),2)
Talla<-`names<-`(data.frame(matrix(c("Height (mts)",paste(Talla.media,"(","±",Talla.sd,")")),ncol = 2,byrow = T)),columns)


#IMC
num1<-c(paste(round(median(base1$imc,na.rm = T ),2),
              paste0('(',round(quantile(base1$imc,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$imc,na.rm = T,probs = c(0.75)),2),')')))

IMC<-`names<-`(data.frame(matrix(c("Body Mass Index (kg/m^2)",num1),ncol = 2)),columns)

#Categorias IMC
IMC.CAT <- table(base2$BMI_CAT,useNA = "always")[c(1:4)]
IMC.CAT.prop <- round(prop.table(table(base2$BMI_CAT,useNA = "always")),4)[c(1:4)]*100
IMC.CAT.FINAL<-`names<-`(data.frame(c("Underweight (%)","Normalweight (%)","Overweight (%)", "Obesity (%)"),
                                    matrix(c(paste(IMC.CAT,paste0('(',IMC.CAT.prop,')'))),ncol = 1,byrow = F)),
                         columns)

#Peso Seco

num1<-c(paste(round(median(base1$peso_seco,na.rm = T ),2),
              paste0('(',round(quantile(base1$peso_seco,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$peso_seco,na.rm = T,probs = c(0.75)),2),')')))

Peso.Seco<-`names<-`(data.frame(matrix(c("Dry Weight(kg)",num1),ncol = 2)),columns)


#Tiempo TSR
num1<-c(paste(round(median(base1$anos_hdcodificado,na.rm = T ),2),
              paste0('(',round(quantile(base1$anos_hdcodificado,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$anos_hdcodificado,na.rm = T,probs = c(0.75)),2),')')))

Tiempo.TSR<-`names<-`(data.frame(matrix(c("Time since RRT (Years)",num1),ncol = 2)),columns)

#Hemoglobina
num1<-c(paste(round(median(base1$hemoglobina,na.rm = T ),2),
              paste0('(',round(quantile(base1$hemoglobina,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$hemoglobina,na.rm = T,probs = c(0.75)),2),')')))

Hemoglobina<-`names<-`(data.frame(matrix(c("Heamoglobin (gr/dl)",num1),ncol = 2)),columns)


#Calcio
num1<-c(paste(round(median(base1$calcio,na.rm = T ),2),
              paste0('(',round(quantile(base1$calcio,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$calcio,na.rm = T,probs = c(0.75)),2),')')))

Calcio<-`names<-`(data.frame(matrix(c("Calcium (mg/dl)",num1),ncol = 2)),columns)


#Sodio
num1<-c(paste(round(median(base1$sodio_serico,na.rm = T ),2),
              paste0('(',round(quantile(base1$sodio_serico,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$sodio_serico,na.rm = T,probs = c(0.75)),2),')')))

Sodio<-`names<-`(data.frame(matrix(c("Sodium (mg/dl)",num1),ncol = 2)),columns)


#Potasio
num1<-c(paste(round(median(base1$potasio,na.rm = T ),2),
              paste0('(',round(quantile(base1$potasio,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$potasio,na.rm = T,probs = c(0.75)),2),')')))

Potasio<-`names<-`(data.frame(matrix(c("Potasium (mg/dl)",num1),ncol = 2)),columns)

#Fosforo
num1<-c(paste(round(median(base1$fosforo,na.rm = T ),2),
              paste0('(',round(quantile(base1$fosforo,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$fosforo,na.rm = T,probs = c(0.75)),2),')')))

Fosforo<-`names<-`(data.frame(matrix(c("Phosphorus (mg/dl)",num1),ncol = 2)),columns)


#Creatinine
num1<-c(paste(round(median(base1$creatinina_serica,na.rm = T ),2),
              paste0('(',round(quantile(base1$creatinina_serica,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$creatinina_serica,na.rm = T,probs = c(0.75)),2),')')))

Creatinine<-`names<-`(data.frame(matrix(c("Creatinine (mg/dl)",num1),ncol = 2)),columns)


#Urea
num1<-c(paste(round(median(base1$urea,na.rm = T ),2),
              paste0('(',round(quantile(base1$urea,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$urea,na.rm = T,probs = c(0.75)),2),')')))

Urea<-`names<-`(data.frame(matrix(c("Urea (mg/dl)",num1),ncol = 2)),columns)

#Diuresis Residual
num1<-c(paste(round(median(base1$diuresis_residual,na.rm = T ),2),
              paste0('(',round(quantile(base1$diuresis_residual,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$diuresis_residual,na.rm = T,probs = c(0.75)),2),')')))

Diruses.Residual<-`names<-`(data.frame(matrix(c("Residual Diuresis (ml/day)",num1),ncol = 2)),columns)

#Albumina
num1<-c(paste(round(median(base1$albumina_serica,na.rm = T ),2),
              paste0('(',round(quantile(base1$albumina_serica,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$albumina_serica,na.rm = T,probs = c(0.75)),2),')')))

Albumin<-`names<-`(data.frame(matrix(c("Albumin (mg/dl)",num1),ncol = 2)),columns)

#Bilirrubina Total
num1<-c(paste(round(median(base1$bilirrubina_totales,na.rm = T ),2),
              paste0('(',round(quantile(base1$bilirrubina_totales,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$bilirrubina_totales,na.rm = T,probs = c(0.75)),2),')')))

BT<-`names<-`(data.frame(matrix(c("Total Bilirubin (mg/dl)",num1),ncol = 2)),columns)


#Billirubina Directa
num1<-c(paste(round(median(base1$bilurrubina_directa,na.rm = T ),2),
              paste0('(',round(quantile(base1$bilurrubina_directa,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$bilurrubina_directa,na.rm = T,probs = c(0.75)),2),')')))

BD<-`names<-`(data.frame(matrix(c("Direct Bilirubin (mg/dl)",num1),ncol = 2)),columns)

#Bilirrubina Indirecta
num1<-c(paste(round(median(base1$bilurrubina_indirecta,na.rm = T ),2),
              paste0('(',round(quantile(base1$bilurrubina_indirecta,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$bilurrubina_indirecta,na.rm = T,probs = c(0.75)),2),')')))

BI<-`names<-`(data.frame(matrix(c("Indirect Bilirubin (mg/dl)",num1),ncol = 2)),columns)

#Fosfatasa Alcalina
num1<-c(paste(round(median(base1$fosfatasa_alcalina,na.rm = T ),2),
              paste0('(',round(quantile(base1$fosfatasa_alcalina,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$fosfatasa_alcalina,na.rm = T,probs = c(0.75)),2),')')))

FA<-`names<-`(data.frame(matrix(c("Phosphatase Alcaline (mg/dl)",num1),ncol = 2)),columns)

#Lactic Deshydrogenase
num1<-c(paste(round(median(base1$desihidrogenasa_lactic,na.rm = T ),2),
              paste0('(',round(quantile(base1$desihidrogenasa_lactic,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base1$desihidrogenasa_lactic,na.rm = T,probs = c(0.75)),2),')')))

DL<-`names<-`(data.frame(matrix(c("Lactic Dehydrogenase (mg/dl)",num1),ncol = 2)),columns)


#Sociodemographic Characteristics
Table3.0<-rbind(Sexo,Edad,Hemodial,Tiempo.TSR,PESO,Talla,IMC,IMC.CAT.FINAL,Diabetes,HAS,CVD,Anemia,Retinopatia,Tabaquismo,Alcoholismo,Neurologic)

#Biochemical and Medication Characteristics
Table4.0<-rbind(EPO,Calcioantago,Betabloq,Tiazidas,Diuretico.Asa,IECA,Alfa.Bloqueador,ARAII.bloqueador,Hierro.Sup,Hipoglucemiantes,Multivitamins,Hemoglobina,Calcio,Sodio,Potasio,Fosforo,Creatinine,Urea,Albumin,BT,BD,BI,FA,DL)

#Rennal Faillure Cause
Table5.0<-rbind(Glomerulopatias,Nefro.HAS,Lupus,Hipoplasia.renal,Diabetes.Nefro,Riñon.Poli,Uropatia.Obs,Nefropatia.Uratos,Nefro.desconocida)


#Stratified by any type of viral infection

columns <- c('Parameter',"Withouth Any-Coinfection (n=1,265)","Any-Viral Coinfection (n=39)")

## Male Female

sexo <- table(base2$genero,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
sexoprop <- round(prop.table(table(base2$genero,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Sexo<-`names<-`(data.frame("Men (%)",
                           matrix(c(paste(sexo,paste0('(',sexoprop,')'))),ncol = 2)),
                columns)
stats::chisq.test(table(base2$genero,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Edad

nortest::ad.test(base2$edad)
num1<-c(paste(round(mean(base2[base2$INFECCION_VIRAL==0,]$edad,na.rm = T),2),
              paste0('(',"±",round(sd(base2[base2$INFECCION_VIRAL==0,]$edad,na.rm = T),2),")")))

num2<-c(paste(round(mean(base2[base2$INFECCION_VIRAL==1,]$edad,na.rm = T),2),
              paste0('(',"±",round(sd(base2[base2$INFECCION_VIRAL==1,]$edad,na.rm = T),2),")")))

Edad<-`names<-`(data.frame(matrix(c("Age (Years)",num1,num2),ncol = 3)),columns)

t.test(edad ~ INFECCION_VIRAL, data = base2)

#Modalidad TSR

hemodialisis <- table(base2$modalidad_REC,base2$INFECCION_VIRAL,useNA = "always")[c(1,2,4,5)]
hemodialisis.prop <- round(prop.table(table(base2$modalidad_REC,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(1,2,4,5)]*100
Hemodial<-`names<-`(data.frame(c("Peritoneal Dialysis (%)","Hemodialysis (%)"),
                               matrix(c(paste(hemodialisis,paste0('(',hemodialisis.prop,')'))),ncol = 2,byrow = F)),
                    columns)

stats::chisq.test(table(base2$modalidad_REC==2,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Glomerulopatias

glomeru <- table(base2$glomerulopatias,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
glomeru.prop <- round(prop.table(table(base2$glomerulopatias,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Glomerulopatias<-`names<-`(data.frame("Glomerulopathies (%)",
                                      matrix(c(paste(glomeru,paste0('(',glomeru.prop,')'))),ncol = 2)),
                           columns)
stats::chisq.test(table(base2$glomerulopatias,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Nefropatia Hipertensiva

nefrop.has <- table(base2$nefropatia_hipertensiva,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
nefrop.has.prop <- round(prop.table(table(base2$nefropatia_hipertensiva,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Nefro.HAS<-`names<-`(data.frame("Hypertensive Nefropathy (%)",
                                matrix(c(paste(nefrop.has,paste0('(',nefrop.has.prop,')'))),ncol = 2)),
                     columns)

stats::chisq.test(table(base2$nefropatia_hipertensiva,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Lupus

lupus <- table(base2$lupus_eritematoso_sistemico,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
lupus.prop <- round(prop.table(table(base2$lupus_eritematoso_sistemico,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Lupus<-`names<-`(data.frame("Lupus (%)",
                            matrix(c(paste(nefrop.has,paste0('(',nefrop.has.prop,')'))),ncol = 2)),
                 columns)

stats::chisq.test(table(base2$lupus_eritematoso_sistemico,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Hypoplastic Nepropathy

hipoplasia <- table(base2$hipoplasia_renal,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
hipoplasia.prop <- round(prop.table(table(base2$hipoplasia_renal,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Hipoplasia.renal<-`names<-`(data.frame("Hypoplastic Nepropathy (%)",
                                       matrix(c(paste(hipoplasia,paste0('(',hipoplasia.prop,')'))),ncol = 2)),
                            columns)

stats::chisq.test(table(base2$hipoplasia_renal,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Diabetes Nefropatia

diab <- table(base2$diabetes_melitus,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
diab.prop <- round(prop.table(table(base2$diabetes_melitus,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Diabetes.Nefro<-`names<-`(data.frame("Diabetic Nepropathy (%)",
                                     matrix(c(paste(diab,paste0('(',diab.prop,')'))),ncol = 2)),
                          columns)

stats::chisq.test(table(base2$diabetes_melitus,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Riñon Poliquistico
poliqu <- table(base2$rinon_poliquistico,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
poliqu.prop <- round(prop.table(table(base2$rinon_poliquistico,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Riñon.Poli<-`names<-`(data.frame("Polycystic Kidney Disease (%)",
                                 matrix(c(paste(poliqu,paste0('(',poliqu.prop,')'))),ncol = 2)),
                      columns)

stats::chisq.test(table(base2$rinon_poliquistico,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Uropatia Obstructiva
obstructiva <- table(base2$uropatia_obstructuva_reflujo,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
obstructiva.prop <- round(prop.table(table(base2$uropatia_obstructuva_reflujo,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Uropatia.Obs<-`names<-`(data.frame("Obstructive Uropathy (%)",
                                   matrix(c(paste(obstructiva,paste0('(',obstructiva.prop,')'))),ncol = 2)),
                        columns)

stats::chisq.test(table(base2$uropatia_obstructuva_reflujo,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Nefropatia Uratos
nefro.uro <- table(base2$nefropatia_por_uratos,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
nefro.uro.prop <- round(prop.table(table(base2$nefropatia_por_uratos,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Nefropatia.Uratos<-`names<-`(data.frame("Urate Nephropathy (%)",
                                        matrix(c(paste(nefro.uro,paste0('(',nefro.uro.prop,')'))),ncol = 2)),
                             columns)

stats::chisq.test(table(base2$nefropatia_por_uratos,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Desconocida
desconocida <- table(base2$desconocido,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
desconocida.prop <- round(prop.table(table(base2$desconocido,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Nefro.desconocida<-`names<-`(data.frame("Unclassified Nephropathy (%)",
                                        matrix(c(paste(desconocida,paste0('(',desconocida.prop,')'))),ncol = 2)),
                             columns)
stats::chisq.test(table(base2$desconocido,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Tabaquismo
tabaq <- table(base2$tabaquismo_actual,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
tabaq.prop <- round(prop.table(table(base2$tabaquismo_actual,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Tabaquismo<-`names<-`(data.frame("Current Smoking (%)",
                                 matrix(c(paste(tabaq,paste0('(',tabaq.prop,')'))),ncol = 2)),
                      columns)

stats::chisq.test(table(base2$tabaquismo_actual,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Alcoholismo
alcoholismo <- table(base2$alcoholismo,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
alcoholismo.prop <- round(prop.table(table(base2$alcoholismo,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Alcoholismo<-`names<-`(data.frame("Alcoholism (%)",
                                  matrix(c(paste(alcoholismo,paste0('(',alcoholismo.prop,')'))),ncol = 2)),
                       columns)

stats::chisq.test(table(base2$alcoholismo,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Diabetes
diabetes <- table(base2$dm,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
diabetes.prop <- round(prop.table(table(base2$dm,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Diabetes<-`names<-`(data.frame("Diabetes (%)",
                               matrix(c(paste(diabetes,paste0('(',diabetes.prop,')'))),ncol = 2)),
                    columns)

stats::chisq.test(table(base2$dm,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Arterial Hypertension
hiper.artr <- table(base2$hta,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
hiper.artr.prop <- round(prop.table(table(base2$hta,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
HAS<-`names<-`(data.frame("Arterial Hypertension (%)",
                          matrix(c(paste(hiper.artr,paste0('(',hiper.artr.prop,')'))),ncol = 2)),
               columns)

stats::chisq.test(table(base2$hta,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Previous Cardiovascular Disease
cardivas <- table(base2$enf_cardiovascular,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
cardivas.prop <- round(prop.table(table(base2$enf_cardiovascular,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
CVD<-`names<-`(data.frame("Previous CVD (%)",
                          matrix(c(paste(cardivas,paste0('(',cardivas.prop,')'))),ncol = 2)),
               columns)

stats::chisq.test(table(base2$enf_cardiovascular,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Retinopatia
retino <- table(base2$retinopatia,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
retino.prop <- round(prop.table(table(base2$retinopatia,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Retinopatia<-`names<-`(data.frame("Retinopathy (%)",
                                  matrix(c(paste(retino,paste0('(',retino.prop,')'))),ncol = 2)),
                       columns)

stats::chisq.test(table(base2$retinopatia,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Alteracion Neurologica
neurologica <- table(base2$alteracion_neurologica,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
neurologica.prop <- round(prop.table(table(base2$alteracion_neurologica,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Neurologic<-`names<-`(data.frame("Neurological Impairment (%)",
                                 matrix(c(paste(retino,paste0('(',retino.prop,')'))),ncol = 2)),
                      columns)

stats::chisq.test(table(base2$alteracion_neurologica,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Anemia
anemia <- table(base2$anemia,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
anemia.prop <- round(prop.table(table(base2$anemia,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Anemia<-`names<-`(data.frame("Anemia (%)",
                             matrix(c(paste(anemia,paste0('(',anemia.prop,')'))),ncol = 2)),
                  columns)

stats::chisq.test(table(base2$anemia,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Eritropoyetina
eritro <- table(base2$utiliza_eritropoyetina,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
eritro.prop <- round(prop.table(table(base2$utiliza_eritropoyetina,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
EPO<-`names<-`(data.frame("Erythropoietin use (%)",
                          matrix(c(paste(eritro,paste0('(',eritro.prop,')'))),ncol = 2)),
               columns)

stats::chisq.test(table(base2$utiliza_eritropoyetina,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Calcioantagonistas
calcio.bloq <- table(base2$calcioantagonista,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
calcio.bloq.prop <- round(prop.table(table(base2$calcioantagonista,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Calcioantago<-`names<-`(data.frame("CCBs use (%)",
                                   matrix(c(paste(calcio.bloq,paste0('(',calcio.bloq.prop,')'))),ncol = 2)),
                        columns)

stats::chisq.test(table(base2$calcioantagonista,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Betabloqueadores
beta.bloq <- table(base2$betabloqueador,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
beta.bloq.prop <- round(prop.table(table(base2$betabloqueador,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Betabloq<-`names<-`(data.frame("Beta-blockers use (%)",
                               matrix(c(paste(beta.bloq,paste0('(',beta.bloq.prop,')'))),ncol = 2)),
                    columns)

stats::chisq.test(table(base2$betabloqueador,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Tiazide Use
tiazida.bloq <- table(base2$diuretico_tiazida,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
tiazida.bloq.prop <- round(prop.table(table(base2$diuretico_tiazida,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Tiazidas<-`names<-`(data.frame("Thiazides use (%)",
                               matrix(c(paste(tiazida.bloq,paste0('(',tiazida.bloq.prop,')'))),ncol = 2)),
                    columns)

stats::chisq.test(table(base2$diuretico_tiazida,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Diureticos de Asa
diuretico.asa <- table(base2$diuretico_de_asa,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
diuretico.asa.prop <- round(prop.table(table(base2$diuretico_de_asa,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Diuretico.Asa<-`names<-`(data.frame("Loop Diuretics use (%)",
                                    matrix(c(paste(diuretico.asa,paste0('(',diuretico.asa.prop,')'))),ncol = 2)),
                         columns)

stats::chisq.test(table(base2$diuretico_de_asa,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#IECA
eca.bloq <- table(base2$eca,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
eca.bloq.prop <- round(prop.table(table(base2$eca,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
IECA<-`names<-`(data.frame("ACEI use (%)",
                           matrix(c(paste(eca.bloq,paste0('(',eca.bloq.prop,')'))),ncol = 2)),
                columns)

stats::chisq.test(table(base2$eca,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Alfa Bloq
alfa.bloq <- table(base2$alfabloqueador,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
alfa.bloq.prop <- round(prop.table(table(base2$alfabloqueador,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Alfa.Bloqueador<-`names<-`(data.frame("Alpha Blockers use (%)",
                                      matrix(c(paste(alfa.bloq,paste0('(',alfa.bloq.prop,')'))),ncol = 2)),
                           columns)

stats::chisq.test(table(base2$alfabloqueador,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#ARA-II Bloq
ARAII.bloq <- table(base2$ara_ii,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
ARAII.prop <- round(prop.table(table(base2$ara_ii,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
ARAII.bloqueador<-`names<-`(data.frame("ARBs use (%)",
                                       matrix(c(paste(ARAII.bloq,paste0('(',ARAII.prop,')'))),ncol = 2)),
                            columns)

stats::chisq.test(table(base2$ara_ii,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
#Hierro Suplementario
hierro <- table(base2$hierro,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
hierro.prop <- round(prop.table(table(base2$hierro,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Hierro.Sup<-`names<-`(data.frame("Iron Supplement Intake (%)",
                                 matrix(c(paste(hierro,paste0('(',hierro.prop,')'))),ncol = 2)),
                      columns)

stats::chisq.test(table(base2$hierro,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Uso Hipoglucemiantes
hipoglu <- table(base2$hipoglucemiantes,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
hipoglu.prop <- round(prop.table(table(base2$hipoglucemiantes,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Hipoglucemiantes<-`names<-`(data.frame("Hypoglycemiant Intake (%)",
                                       matrix(c(paste(hipoglu,paste0('(',hipoglu.prop,')'))),ncol = 2)),
                            columns)

stats::chisq.test(table(base2$hipoglucemiantes,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Multivitamins
multivit <- table(base2$multivitaminas,base2$INFECCION_VIRAL,useNA = "always")[c(2,5)]
multivit.prop <- round(prop.table(table(base2$multivitaminas,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
Multivitamins<-`names<-`(data.frame("Hypoglycemiant Intake (%)",
                                    matrix(c(paste(multivit,paste0('(',multivit.prop,')'))),ncol = 2)),
                         columns)

stats::chisq.test(table(base2$multivitaminas,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Peso
num1<-c(paste(round(mean(base2[base2$INFECCION_VIRAL==0,]$peso_actual,na.rm = T),2),
              paste0('(',"±",round(sd(base2[base2$INFECCION_VIRAL==0,]$peso_actual,na.rm = T),2),")")))

num2<-c(paste(round(mean(base2[base2$INFECCION_VIRAL==1,]$peso_actual,na.rm = T),2),
              paste0('(',"±",round(sd(base2[base2$INFECCION_VIRAL==1,]$peso_actual,na.rm = T),2),")")))

PESO<-`names<-`(data.frame(matrix(c("Weight (kg)",num1,num2),ncol = 3)),columns)

wilcox.test(peso_actual~INFECCION_VIRAL,data = base2)

#Talla
num1<-c(paste(round(mean(base2[base2$INFECCION_VIRAL==0,]$talla,na.rm = T),2),
              paste0('(',"±",round(sd(base2[base2$INFECCION_VIRAL==0,]$talla,na.rm = T),2),")")))

num2<-c(paste(round(mean(base2[base2$INFECCION_VIRAL==1,]$talla,na.rm = T),2),
              paste0('(',"±",round(sd(base2[base2$INFECCION_VIRAL==1,]$talla,na.rm = T),2),")")))

Talla<-`names<-`(data.frame(matrix(c("Height (mts)",num1,num2),ncol = 3)),columns)

#IMC
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$imc,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$imc,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$imc,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$imc,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$imc,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$imc,na.rm = T,probs = c(0.75)),2),')')))

IMC<-`names<-`(data.frame(matrix(c("Body Mass Index (kg/m^2)",num1,num2),ncol = 3)),columns)

wilcox.test(imc~INFECCION_VIRAL,data = base2)
#Categorias IMC
IMC.CAT <- table(base2$BMI_CAT,base2$INFECCION_VIRAL,useNA = "always")[c(1,2,3,4,6,7,8,9)]
IMC.CAT.prop <- round(prop.table(table(base2$BMI_CAT,base2$INFECCION_VIRAL,useNA = "always"),2),4)[c(2,5)]*100
IMC.CAT.FINAL<-`names<-`(data.frame(c("Underweight (%)","Normalweight (%)","Overweight (%)", "Obesity (%)"),
                                    matrix(c(paste(IMC.CAT,paste0('(',IMC.CAT.prop,')'))),ncol = 2)),
                         columns)

stats::chisq.test(table(base2$BMI_CAT==1,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
stats::chisq.test(table(base2$BMI_CAT==2,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
stats::chisq.test(table(base2$BMI_CAT==3,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)
stats::chisq.test(table(base2$BMI_CAT==4,base2$INFECCION_VIRAL),correct = T,simulate.p.value = T)

#Peso Seco
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$peso_seco,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$peso_seco,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$peso_seco,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$peso_seco,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$peso_seco,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$peso_seco,na.rm = T,probs = c(0.75)),2),')')))

Peso.Seco<-`names<-`(data.frame(matrix(c("Dry Weight(kg)",num1,num2),ncol = 3)),columns)

wilcox.test(peso_seco~INFECCION_VIRAL,data = base2)

#Tiempo TSR
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$anos_hdcodificado,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$anos_hdcodificado,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$anos_hdcodificado,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$anos_hdcodificado,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$anos_hdcodificado,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$anos_hdcodificado,na.rm = T,probs = c(0.75)),2),')')))

Tiempo.TSR<-`names<-`(data.frame(matrix(c("Time since RRT (Years)",num1,num2),ncol = 3)),columns)

wilcox.test(anos_hdcodificado~INFECCION_VIRAL,data = base2)
#Hemoglobina
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$hemoglobina,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$hemoglobina,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$hemoglobina,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$hemoglobina,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$hemoglobina,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$hemoglobina,na.rm = T,probs = c(0.75)),2),')')))

Hemoglobina<-`names<-`(data.frame(matrix(c("Heamoglobin (gr/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(hemoglobina~INFECCION_VIRAL,data = base2)

#Calcio
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$calcio,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$calcio,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$calcio,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$calcio,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$calcio,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$calcio,na.rm = T,probs = c(0.75)),2),')')))

Calcio<-`names<-`(data.frame(matrix(c("Calcium (gr/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(calcio~INFECCION_VIRAL,data = base2)

#Sodio
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$sodio_serico,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$sodio_serico,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$sodio_serico,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$sodio_serico,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$sodio_serico,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$sodio_serico,na.rm = T,probs = c(0.75)),2),')')))

Sodio<-`names<-`(data.frame(matrix(c("Sodium (gr/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(sodio_serico~INFECCION_VIRAL,data = base2)

#Potasio
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$potasio,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$potasio,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$potasio,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$potasio,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$potasio,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$potasio,na.rm = T,probs = c(0.75)),2),')')))

Potasio<-`names<-`(data.frame(matrix(c("Potasium (gr/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(potasio~INFECCION_VIRAL,data = base2)
#Fosforo

num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$fosforo,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$fosforo,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$fosforo,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$fosforo,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$fosforo,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$fosforo,na.rm = T,probs = c(0.75)),2),')')))

Fosforo<-`names<-`(data.frame(matrix(c("Phosphorus (gr/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(fosforo~INFECCION_VIRAL,data = base2)
#Creatinine
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$creatinina_serica,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$creatinina_serica,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$creatinina_serica,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$creatinina_serica,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$creatinina_serica,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$creatinina_serica,na.rm = T,probs = c(0.75)),2),')')))

Creatinine<-`names<-`(data.frame(matrix(c("Creatinine (gr/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(creatinina_serica~INFECCION_VIRAL,data = base2)
#Urea
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$urea,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$urea,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$urea,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$urea,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$urea,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$urea,na.rm = T,probs = c(0.75)),2),')')))

Urea<-`names<-`(data.frame(matrix(c("Urea (gr/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(urea~INFECCION_VIRAL,data = base2)
#Diuresis Residual
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$diuresis_residual,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$diuresis_residual,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$diuresis_residual,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$diuresis_residual,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$diuresis_residual,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$diuresis_residual,na.rm = T,probs = c(0.75)),2),')')))

Diruses.Residual<-`names<-`(data.frame(matrix(c("Residual Diuresis (ml/day)",num1,num2),ncol = 3)),columns)

wilcox.test(diuresis_residual~INFECCION_VIRAL,data = base2)
#Albumina
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$albumina_serica,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$albumina_serica,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$albumina_serica,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$albumina_serica,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$albumina_serica,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$albumina_serica,na.rm = T,probs = c(0.75)),2),')')))

Albumin<-`names<-`(data.frame(matrix(c("Albumin (mg/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(albumina_serica~INFECCION_VIRAL,data = base2)
#Bilirrubina Total
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$bilirrubina_totales,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$bilirrubina_totales,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$bilirrubina_totales,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$bilirrubina_totales,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$bilirrubina_totales,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$bilirrubina_totales,na.rm = T,probs = c(0.75)),2),')')))

BT<-`names<-`(data.frame(matrix(c("Total Bilirubin (mg/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(bilirrubina_totales~INFECCION_VIRAL,data = base2)
#Billirubina Directa
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$bilurrubina_directa,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$bilurrubina_directa,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$bilurrubina_directa,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$bilurrubina_directa,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$bilurrubina_directa,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$bilurrubina_directa,na.rm = T,probs = c(0.75)),2),')')))

BD<-`names<-`(data.frame(matrix(c("Direct Bilirubin (mg/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(bilurrubina_directa~INFECCION_VIRAL,data = base2)

#Bilirrubina Indirecta
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$bilurrubina_indirecta,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$bilurrubina_indirecta,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$bilurrubina_indirecta,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$bilurrubina_indirecta,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$bilurrubina_indirecta,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$bilurrubina_indirecta,na.rm = T,probs = c(0.75)),2),')')))

BI<-`names<-`(data.frame(matrix(c("Indirect Bilirubin (mg/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(bilurrubina_indirecta~INFECCION_VIRAL,data = base2)
#Fosfatasa Alcalina
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$fosfatasa_alcalina,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$fosfatasa_alcalina,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$fosfatasa_alcalina,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$fosfatasa_alcalina,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$fosfatasa_alcalina,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$fosfatasa_alcalina,na.rm = T,probs = c(0.75)),2),')')))

FA<-`names<-`(data.frame(matrix(c("Phosphatase Alcaline (mg/dl)",num1,num2),ncol = 3)),columns)

wilcox.test(fosfatasa_alcalina~INFECCION_VIRAL,data = base2)
#Lactic Deshydrogenase
num1<-c(paste(round(median(base2[base2$INFECCION_VIRAL==0,]$desihidrogenasa_lactic,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==0,]$desihidrogenasa_lactic,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==0,]$desihidrogenasa_lactic,na.rm = T,probs = c(0.75)),2),')')))

num2<-c(paste(round(median(base2[base2$INFECCION_VIRAL==1,]$desihidrogenasa_lactic,na.rm = T ),2),
              paste0('(',round(quantile(base2[base2$INFECCION_VIRAL==1,]$desihidrogenasa_lactic,na.rm = T,probs = c(0.25)),2),"-",
                     round(quantile(base2[base2$INFECCION_VIRAL==1,]$desihidrogenasa_lactic,na.rm = T,probs = c(0.75)),2),')')))

DL<-`names<-`(data.frame(matrix(c("Lactic Dehydrogenase (mg/dl)",num1,num2),ncol = 3)),columns)

#Sociodemographic Characteristics
Table3.1<-rbind(Sexo,Edad,Hemodial,Tiempo.TSR,PESO,Talla,IMC,IMC.CAT.FINAL,Diabetes,HAS,CVD,Anemia,Retinopatia,Tabaquismo,Alcoholismo,Neurologic)
Table3<-cbind(Table3.0,Table3.1[,c(-1)])
Table3_Flex<-flextable::align(flextable::flextable(Table3,cwidth=4),align="center",part="all")%>%flextable::autofit()
#flextable::save_as_docx(Table3_Flex,path="Table_1.docx")

#Biochemical and Medication Characteristics
Table4.1<-rbind(EPO,Calcioantago,Betabloq,Tiazidas,Diuretico.Asa,IECA,Alfa.Bloqueador,ARAII.bloqueador,Hierro.Sup,Hipoglucemiantes,Multivitamins,Hemoglobina,Calcio,Sodio,Potasio,Fosforo,Creatinine,Urea,Albumin,BT,BD,BI,FA,DL)
Table4<-cbind(Table4.0,Table4.1[,c(-1)])
Table4_Flex<-flextable::align(flextable::flextable(Table4,cwidth=4),align="center",part="all")%>%flextable::autofit()
#flextable::save_as_docx(Table4_Flex,path="Table_2.docx")

#Rennal Faillure Cause
Table5.1<-rbind(Glomerulopatias,Nefro.HAS,Lupus,Hipoplasia.renal,Diabetes.Nefro,Riñon.Poli,Uropatia.Obs,Nefropatia.Uratos,Nefro.desconocida)
Table5<-cbind(Table5.0,Table5.1[,c(-1)])
Table5_Flex<-flextable::align(flextable::flextable(Table5,cwidth=4),align="center",part="all")%>%flextable::autofit()
#flextable::save_as_docx(Table5_Flex,path="Table_3.docx")







