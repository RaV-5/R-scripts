library(tidyverse)
library(readxl)

#loading .xls file + /dir ??
filename <- excel_sheets("test.xls")

fun <- function(filename){
  for (x in filename){
    if (x == "Results"){
      plik <- read_excel("test.xls", sheet="Results", skip=42)
      break
    }
    else if ((x==tail(excel_sheets("test.xls"), n=1))){
      print("Brak zakladki 'Results' w pliku excel")
    }
  }
  
  #global variables
  # ref_gene <- readline(prompt="Podaj gen referencyjny: ")
  ref_gene <- "B-actin"
  
  #data cleaning
  data <- plik %>%
    select("name"="Sample Name", "target"="Target Name", "ct"="CT")%>%
    filter(name != "NTC" & ct!="Undetermined") %>%
    mutate(ct=as.numeric(ct)) %>%
    drop_na()
  
  #calculating mean of ref_gene
  kal_ref <- data %>%
    filter(target==ref_gene)%>%
    group_by(name, target)%>%
    summarise(ct_ref=mean(ct))%>%
    ungroup()%>%
    select(!target)
  
  #joining ref_gene CT to all samples according to names and calculating delta ct
  delta <- data %>%
    right_join(y=kal_ref, by=c("name"))%>%
    mutate("delta_ct"= ct-ct_ref)%>%
    filter(target!=ref_gene)%>%
    select(name, target, delta_ct)
  
  #extracting only calibrator from data
  cal <- delta %>%
    filter(name=="kalibrator")%>%
    select(target, delta_ct)%>% 
    rename(delta_kal=delta_ct)%>%
    group_by(target)%>%
    summarise(delta_kal=mean(delta_kal))
  
  #joining calibrator to table and calculating delta delta ct 
  dd <- delta %>%
    left_join(y=cal, by=c("target"))%>%
    mutate(dd=delta_ct-delta_kal)%>%
    filter(name!="kalibrator")%>%
    mutate(dd2 = log10(2^(-dd)))%>%
    group_by(name, target)%>%
    mutate(stdev=sd(dd2), dd2_mean=mean(dd2))
  
  ## zapisanie wykresu --> mozna wyrzucic paste0 i zapisze automatycznie w wd
  wd=as.character(getwd())
  bmp(file= paste0(wd, "/plot.bmp"))
  
  #plotting
  plot <- ggplot(dd, aes(x=name, y=dd2_mean, fill=target)) +
    geom_col(position="dodge")+
    labs(y= "log(dd2)", x="Sample Name", fill="Target gene")+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "White"))+
    scale_fill_grey()+
    geom_hline(yintercept=0)+
    geom_errorbar(aes(ymin=dd2_mean-stdev, 
                      ymax=dd2_mean +stdev), 
                  position="dodge")
  print(plot)
  dev.off()
}

fun(filename)

