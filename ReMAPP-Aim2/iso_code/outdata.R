#*****************************************************************************
#Results data function
#Author: Xiaoyan Hu
#Email; xyh@gwu.edu
#*****************************************************************************
outdata <- function(data, outcome) {
df <- as.data.frame(data$groups) %>% 
  bind_cols(as.data.frame(data$estimates)) %>% 
  group_by(`data$groups`) %>%  
  summarise(
    group_n = n(),
    risk_mscore = sum(`data$estimates`)/n())

df_out <- as.data.frame(
  list(
    Outcome = outcome,
    N1 = df$group_n[1],
    "Risk/Mscore1" = df$risk_mscore[1],
    Thres1 = data$brkPoints[1],
    N2 = df$group_n[2],
    "Risk/Mscore2" = df$risk_mscore[2],
    Thres2 = data$brkPoints[2],
    N3 = df$group_n[3],
    "Risk/Mscore3" = df$risk_mscore[3],
    Thres3 = data$brkPoints[3],
    N4 = df$group_n[4],
    "Risk/Mscore4" = df$risk_mscore[4],
    Thres4 = data$brkPoints[4],
    N5 = df$group_n[5],
    "Risk/Mscore5" = df$risk_mscore[5],
    Thres5 = data$brkPoints[5],
    N6 = df$group_n[6],
    "Risk/Mscore6" = df$risk_mscore[6],
    Thres6 = data$brkPoints[6],
    N7 = df$group_n[7],
    "Risk/Mscore7" = df$risk_mscore[7],
    Thres7 = data$brkPoints[7],
    N8 = df$group_n[8],
    "Risk/Mscore8" = df$risk_mscore[8],
    Thres8 = data$brkPoints[8],
    N9 = df$group_n[9],
    "Risk/Mscore9" = df$risk_mscore[9]
  ))
}



