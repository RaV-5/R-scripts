library('tidyverse')
library('readxl')
library('ggpubr')

#data <- read_excel('2022-11-18 HL-60_Kasumi-1_GAPDH.xlsx', skip = 42, sheet = 'Results')
#data <- data[-20,]
#data <- read_excel('2022-10-28 Jurkat_MV_GAPDH_BA2.xlsx', skip = 42, sheet = 'Results')
#data <- read_excel('2022-10-31_HEK_GAPDH.xlsx', skip = 40, sheet = 'Results')
#data <- read_excel('2022-12-15_electroporation HL6- Jurkat.xls', skip = 42, sheet = 'Results')
#data <- read_excel('2022-12-19_reversed transfection_NB4_Kasumi-1.xlsx', skip = 40, sheet = 'Results')
#data <- read_excel('2022-12-29_Kasumi-1_electroporation_GAPDH.xlsx', skip = 40, sheet = 'Results')
data <- read_excel('Jurkat.xlsx', skip = 42, sheet = 'Results')

ref_gen <- 'B-actin'
method <- 'electroporation'
gene <- 'miR-193b'


data <- data %>%
  separate('Sample Name', c('Cell_line', 'Vargroup'), '_')

for (x in levels(as.factor(data$Cell_line))){
  if (x != 'K+' && x!= 'NTC'){
    
    
# DODAĆ AUTOMATYCZNE WYKLUCZANIE ODSTAJACYCH POMIAROW
    
    
    #data preprocessing + mean and SD for technical repeats, Undetermined CT are discarded
    cells <- data %>%
      filter(Cell_line == x) %>%
      separate('Well Position', c('Plate_row', 'Plate_col'), 1) %>%
      select(Plate_col, Plate_row, Cell_line, Vargroup, Target = 'Target Name', CT = 'CT') %>%
      filter(CT != 'Undetermined') %>%
      mutate(CT = as.double(CT)) %>%
      group_by(Vargroup, Target, Plate_col) %>%
      summarise(Cell_line, CT, Plate_row, mean_CT = mean(CT), SD = sd(CT, na.rm = TRUE))
    
    Exp_groups <- unique(pull(cells, Vargroup))
    
#filter for ref gene and table join for delta CT calculations [gene - ref]
    ref <- cells %>%
      filter(Target == ref_gen) %>%
      distinct(Plate_col, Vargroup, Target, mean_CT, SD)

    dCT <- cells %>%
      filter(Target != ref_gen) %>%
      right_join(y = ref, by = c('Vargroup', 'Plate_col'), suffix = c('', '_ref')) %>%
      mutate(dCT = mean_CT - mean_CT_ref) %>%
      select(-c(CT, Plate_row)) %>%
      group_by(Vargroup, Plate_col) %>%
      mutate(mean_dCT = mean(dCT, na.rm = TRUE))%>%
      distinct() %>%
      na.omit()      # dodano OMIT, zwracać uwagę czy nie wpływa na wyniki!!!

#ddCT
    refCT <- mean(dCT[dCT$Vargroup == 'K-', "mean_dCT"][[1]])
    ddCT <- dCT %>%
      mutate(ref_sample = refCT) %>%
      mutate(ddCT = mean_dCT - ref_sample, final = 2^(-(ddCT)))

#plotting preprocessing + levels w factor zmienia kolejnosc kolumn na wykresie
    plot <- ddCT %>%
      group_by(Vargroup) %>%
      mutate(Vargroup = factor(Vargroup, levels= Exp_groups)) %>%
      mutate(Cell_line, SD = mean(SD, na.rm = T), final = mean(final)) %>%
      distinct(Vargroup, Cell_line, SD, final)

#statistic test t-test WILCOX VS T-TEST
    com <- compare_means(data = ddCT, method = 't.test', formula = dCT~Vargroup, ref.group = 'K-')
    com

    exp_gene <- unique(ddCT$Target)

#plotting 
    xx <- ggplot(plot, aes(y = final, x = Vargroup))+
      geom_col(position = 'dodge', 
               color = 'black', 
               fill = 'grey', 
               width = 0.5)+
      labs(y= "Relative expression (log scale)", 
           x="Sample Name", 
           fill="Target gene")+
      geom_errorbar(aes(ymin=final - SD, 
                        ymax=final + SD), 
                    width = 0.1)+
      scale_y_continuous(trans = "log10")+
      expand_limits(y = c(0.1, 3))+
      geom_hline(yintercept=1, size = 0.8)+
      theme_classic()+
      theme(panel.border = element_rect(colour = "black", 
                                       fill=NA, size=1.4),
      plot.title = element_text(hjust = 0.5))+
      labs(title = sprintf('Relative expression of %s in %s cell line.', 
                           exp_gene, x)) +
      stat_pvalue_manual(data = com, y.position = c(0.25, 0.45), 
                         label = 'T-test, p = {p.format}')
    
    print(xx)

    ggsave(sprintf('Plot_%s_%s_%s.jpg', method, x, Sys.Date()), xx, width = 10, height = 10, dpi = 400)

  }
  else {print('Not a cell line.[K+/NTC]')}
}

