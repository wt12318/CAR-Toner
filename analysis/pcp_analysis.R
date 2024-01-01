library(dplyr)
vhh <- data.table::fread("data/VHH patent_sequence",data.table = F)
write.csv(vhh %>% select(sequence) %>% rename(seq=sequence),
          file = "data/vhh_pre.csv",row.names = F,quote = F)
vhh_res <- read.csv("data/vhh_pre_res.csv") %>% select(-X)
vhh_res$pre <- round(2^(vhh_res$pre))

VLR <- read.csv("data/VLR_sequence.csv",header = T)
VLR2 <- readxl::read_xlsx("data/uniprot-VLR sequence.xlsx")
VLR <- bind_rows(VLR %>% select(Entry.ID,Sequence) %>% rename(Entry=Entry.ID),
                 VLR2 %>% select(Entry,Sequence)) %>% distinct_all() %>% 
  filter(nchar(Sequence)>1)
VLR <- VLR[!duplicated(VLR$Sequence),]
write.csv(VLR %>% rename(seq=Sequence),
          file = "data/vlr_pre.csv",row.names = F,quote = F)
vlr_res <- read.csv("data/vlr_pre_res.csv") %>% select(-X)
vlr_res$pre <- round(2^(vlr_res$pre))
vlr_res <- vlr_res %>% 
  filter(nchar(seq)>100 & nchar(seq)<300)

shark <- readxl::read_xlsx("data/shark IgNAR.xlsx") %>% select(Sequence)
shark <- shark[!duplicated(shark$Sequence),]
write.csv(shark %>% rename(seq=Sequence),
          file = "data/shark_pre.csv",row.names = F,quote = F)
shark_res <- read.csv("data/shark_pre_res.csv") %>% select(-X)
shark_res$pre <- round(2^(shark_res$pre))

###翻译 scFv
library(seqinr)
get_cdrfr <- function(seq){
  aa_seq <- seqinr::translate(seq) %>% paste(.,collapse = "")
  temp_file <- paste0("/home/wt/PCP_model/tempdir/test.fa")
  seqinr::write.fasta(aa_seq,names = "test",file.out = temp_file)
  commd <- paste0("/home/wt/PCP_model/AbRSA/AbRSA -i ", temp_file, 
                  " -g -o /home/wt/PCP_model/tempdir/test_out",
                  " > /home/wt/PCP_model/tempdir/out_cdrfr")
  system(commd)
  out <- data.table::fread(paste0("/home/wt/PCP_model/tempdir/out_cdrfr"),
                           skip = "#similarity",
                           data.table = F,sep = ":",header = F)
  if (nrow(out) == 1){
    return(NA)
  }else{
    out <- out[2:(nrow(out)-1),]
    out_seq <- paste(out$V2,collapse = "")
    return(out_seq)
  }
}

scfv <- read.fasta("data/scFv sequence.txt",as.string = TRUE,
                   seqtype = "DNA")
scfv_aa <- data.frame(seq=c(1:470))
for (i in 1:nrow(scfv_aa)){
  scfv_aa$seq[i] <- get_cdrfr(s2c(as.character(scfv[[i]])))
}
scfv_aa <- scfv_aa %>% 
  filter(!is.na(seq)) %>% 
  distinct()
write.csv(scfv_aa,
          file = "data/scfv_aa_pre.csv",row.names = F,quote = F)

scfv_res <- read.csv("data/scfv_aa_pre_res.csv") %>% select(-X)
scfv_res$pre <- round(2^(scfv_res$pre))

cardb_seqs_adduser <- readRDS("data/cardb_seqs_adduser.rds")
write.csv(cardb_seqs_adduser,
          file = "data/scfv_db_pre.csv",row.names = F,quote = F)
scfv_db <- read.csv("data/scfv_db_pre.csv") %>% select(-X)
scfv_db$pre <- round(2^(scfv_db$pre))
scfv_res <- bind_rows(
  scfv_res,
  scfv_db %>% select(seq,pre)
)
scfv_res <- scfv_res %>% distinct_all()

all_pcp <- bind_rows(
  scfv_res %>% mutate(type="scfv"),
  shark_res %>% mutate(type="shark"),
  vhh_res %>% mutate(type="vhh"),
  vlr_res %>% mutate(type="vlr") %>% select(-Entry)
)
saveRDS(all_pcp,file = "data/all4pcp.rds")

###
all_pcp <- readRDS("data/all4pcp.rds")
all_pcp <- all_pcp %>% 
  filter(nchar(seq) >100 & nchar(seq)<300) %>% 
  mutate(type = ifelse(type=="shark","vnar",type))
all_pcp <- all_pcp %>% mutate(len=nchar(seq))
scfv <- all_pcp %>% filter(type=="scfv")
vhh <- all_pcp %>% filter(type=="vhh")
vlr <- all_pcp %>% filter(type=="vlr")
vnar <- all_pcp %>% filter(type=="vnar")

all_pcp <- all_pcp %>% 
  mutate(type = toupper(type)) %>% 
  mutate(type = ifelse(type == "SCFV","scFV",type))
all_pcp$type <- factor(all_pcp$type,levels = c("scFV","VHH","VLR","VNAR"))
all_pcp_summ <- all_pcp %>% 
  group_by(type) %>% 
  summarise(per_80 = paste0("Proportion > 120: ",round(mean(pre > 120),4) * 100, "%"))

library(ggtext)
library(ggstatsplot)
library(ggplot2)
ggbetweenstats(
  data  = all_pcp,
  x     = type,
  y     = pre,
  type  = "nonparametric",
  centrality.type = "parametric",
  xlab = "",
  ylab = "PCP score",
  plot.type = "box",
  results.subtitle = FALSE,
  pairwise.display = "all",
  pairwise.comparisons = FALSE
) +
  ggplot2::scale_color_manual(values = c("#bc3c45", "#35a3df", "#fab330","#c6c6c6"))+
  # geom_text_repel(data  = all_pcp, aes(x = type, y = pre, label = name),
  #                 box.padding = 3, max.overlaps = Inf, seed = 123)+
  guides()+
  theme_classic()+
  theme(legend.position = "none")+
  geom_hline(yintercept=46, linetype="dashed", color = "red")+
  geom_hline(yintercept=56, linetype="dashed", color = "red")+
  annotate(geom = "richtext", x = 4, y = 180, 
           label = "<br><span style = 'color:#bc3c45;'>scFV</span>-<span style = 'color:#35a3df;'>VHH</span>:P<sub><i>adj</i></sub> = 2.11e-09
                    <br><span style = 'color:#bc3c45;'>scFV</span>-<span style = 'color:#fab330;'>VLR</span>:P<sub><i>adj</i></sub> = 2.98e-40
                    <br><span style = 'color:#bc3c45;'>scFV</span>-<span style = 'color:#c6c6c6;'>VNAR</span>:P<sub><i>adj</i></sub> = 3.91e-48
                    <br><span style = 'color:#35a3df;'>VHH</span>-<span style = 'color:#fab330;'>VLR</span>:P<sub><i>adj</i></sub> = 1.62e-70
                    <br><span style = 'color:#35a3df;'>VHH</span>-<span style = 'color:#c6c6c6;'>VNAR</span>:P<sub><i>adj</i></sub> = 6.03e-41
                    <br><span style = 'color:#fab330;'>VLR</span>-<span style = 'color:#c6c6c6;'>VNAR</span>:P<sub><i>adj</i></sub> = 5.73e-14"
  )
ggsave(filename = "figs/pcp_diff4_3.pdf",width = 10,height = 8)

###
###his
p1 <- gghistogram(data=scfv,x="pre",fill = "#bc3c45",xlab = "PCP Score",ylab = "Counts",
                  title = "scFV")+
  geom_vline(xintercept=46, linetype="dashed", color = "red",size=1)+
  geom_vline(xintercept=56, linetype="dashed", color = "red",size=1)+
  annotate(geom = "text",x=51,y=65,
           label=paste0(round(mean(scfv$pre >= 46 & scfv$pre <= 56),3) * 100,"%"))

p2 <- gghistogram(data=vhh,x="pre",fill = "#35a3df",xlab = "PCP Score",ylab = "Counts",
                  title = "VHH")+
  geom_vline(xintercept=46, linetype="dashed", color = "red",size=1)+
  geom_vline(xintercept=56, linetype="dashed", color = "red",size=1)+
  annotate(geom = "text",x=51,y=1800,
           label=paste0(round(mean(vhh$pre >= 46 & vhh$pre <= 56),3) * 100,"%"))

p3 <- gghistogram(data=vlr,x="pre",fill = "#fab330",xlab = "PCP Score",ylab = "Counts",
                  title = "VLR")+
  geom_vline(xintercept=46, linetype="dashed", color = "red",size=1)+
  geom_vline(xintercept=56, linetype="dashed", color = "red",size=1)+
  annotate(geom = "text",x=51,y=550,
           label=paste0(round(mean(vlr$pre >= 46 & vlr$pre <= 56),3) * 100,"%"))

p4 <- gghistogram(data=vnar,x="pre",fill = "#c6c6c6",xlab = "PCP Score",ylab = "Counts",
                  title = "VNAR")+
  geom_vline(xintercept=46, linetype="dashed", color = "red",size=1)+
  geom_vline(xintercept=56, linetype="dashed", color = "red",size=1)+
  annotate(geom = "text",x=51,y=60,
           label=paste0(round(mean(vnar$pre >= 46 & vnar$pre <= 56),3) * 100,"%"))
library(patchwork)
(p1 + p2 ) /(p3 + p4)
ggsave("figs/his_compare_2.pdf",width = 15,height = 10)

###
mut_pre <- read.csv("data/cardb_mut_val_add.csv") %>% select(-X)
mut_pre$pre <- round(2^(mut_pre$pre))
mut_pre <- mut_pre %>% 
  mutate(mut_type = gsub("_.+","",id)) %>% 
  mutate(mut_type = case_when(
    mut_type == "kq" ~ "K-Q",
    mut_type == "qk" ~ "Q-K",
    mut_type == "wt" ~ "WT"
  ))
mut_pre$mut_type <- factor(mut_pre$mut_type,levels = c("K-Q","WT","Q-K"))

p1 <- ggboxplot(mut_pre,x="mut_type", y="pre", fill = "mut_type",
                xlab = "Mutation Type",ylab ="PCP Score",
                ggtheme = theme_bw())+
  rotate_x_text(90)
facet(p1, facet.by = "name", ncol = 9)+
  theme(legend.position = "none")
ggsave("figs/all_car_mut_compare.pdf",width = 18,height = 19)


