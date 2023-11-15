library(dplyr)
library(seqinr)

##pdb sequence
all_fasta <- read.fasta("data/pdb_seqres.txt",as.string = TRUE,
                        seqtype = "AA")
##single chain pdb sequence
all_names <- data.frame(name = names(all_fasta)) %>% 
  mutate(pro = gsub("_.+","",name))
all_names_summ <- all_names %>% 
  group_by(pro) %>% 
  summarise(counts=n()) %>% 
  filter(counts == 1)
all_names <- all_names %>% filter(pro %in% all_names_summ$pro)
mono <- all_fasta[all_names$name]

###length in 100 ~ 400
all_lens <- sapply(mono,function(x){nchar(x)},simplify = TRUE)
need_lens <- all_lens[all_lens>100 & all_lens<400]
need_fasta <- names(need_lens) %>% unique()
need_fasta <- paste0("pdb",need_fasta)

need_pro <- gsub("_.+","",need_fasta) %>% unique()
all_pdb <- list.files("/home/data/sda/wt/pdb/") ##downloaded pdb files
need_pro <- need_pro[which(paste0(need_pro,".ent") %in% all_pdb)]
need_fasta <- need_fasta[which(paste0(need_pro,".ent") %in% all_pdb)]

write.table(need_fasta,file = "/home/data/sda/wt/test/need_fa",sep = "\t",quote = F,
            row.names = F,col.names = F)

###calculate PCP based on PDB 
###pcp res
library(parallel)
files <- list.files("/home/data/sda/wt/pcp_res/") ##calculated results
res <- mclapply(files,
                function(x){
                  dt <- read.csv(paste0("/home/data/sda/wt/pcp_res/",x))
                  dt <- data.frame(name=gsub("_patches.csv","",x),pcp=nrow(dt))
                  return(dt)
                }, mc.cores = 40)
res <- bind_rows(res)
res$name <- substr(res$name,4,9)
all_fasta <- read.fasta("data/pdb_seqres.txt",as.string = TRUE,
                        seqtype = "AA")
all_fasta <- all_fasta[res$name]
res$seq <- as.character(all_fasta)
dup_seq <- res$seq[which(duplicated(res$seq))] ##remove duplicates
res <- res %>% filter(!(seq %in% dup_seq))
res <- res %>% filter(pcp > 10)
saveRDS(res,file = "data/pdb_pcp.rds")

####alphafold
uniport <- seqinr::read.fasta("/home/data/sda/wt/uniprot_sprot.fasta.gz",
                              as.string = TRUE,
                              seqtype = "AA")
names(uniport) <- gsub("sp[|]","",names(uniport)) %>% gsub("[|].+","",.)

all_lens <- sapply(uniport,function(x){nchar(x)},simplify = TRUE)
need_lens <- all_lens[all_lens>100 & all_lens<400]
need_id <- names(need_lens) %>% unique()
ids <- paste0("AF-",need_id,"-F1-model_v4.pdb.gz")
ids_pdb <- list.files("/home/data/sda/wt/alphafold_pdb/") ##downloaded
ids <- ids[which(ids %in% ids_pdb)]
ids <- gsub("-F1-model_v4","",ids) %>% gsub("AF-","",.) %>% gsub(".pdb.gz","",.)
need_seq <- uniport[ids]

alpha_pdb <- data.frame(name=ids,seq=as.character(need_seq))
dup_seq <- alpha_pdb$seq[which(duplicated(alpha_pdb$seq))]
alpha_pdb <- alpha_pdb %>% filter(!(seq %in% dup_seq))

##calculate PCP
files <- list.files("/home/data/sda/wt/pcp_res_alpha/")
res <- mclapply(files,
                function(x){
                  
                  tryCatch(
                    {
                      dt <- read.csv(paste0("/home/data/sda/wt/pcp_res_alpha/",x))
                      dt <- data.frame(name=gsub("_patches.csv","",x),pcp=nrow(dt))
                      return(dt)
                    },
                    error = function(cnd) {
                      return(x)
                    }
                  )
                }, mc.cores = 40)

lengths(res) -> res_len

res <- bind_rows(res)
res$name <- gsub("-F1-model_v4","",res$name) %>% 
  gsub("AF-","",.)

uniport <- seqinr::read.fasta("/home/data/sda/wt/uniprot_sprot.fasta.gz",
                              as.string = TRUE,
                              seqtype = "AA")
names(uniport) <- gsub("sp[|]","",names(uniport)) %>% gsub("[|].+","",.)
need_fasta <- uniport[res$name]

res$seq <- as.character(need_fasta)
dup_seq <- res$seq[which(duplicated(res$seq))]
res <- res %>% filter(!(seq %in% dup_seq))
res <- res %>% filter(pcp > 10)
saveRDS(res,file = "data/uniport_pcp.rds")

##merge
uniport <- readRDS("~/PCP_model/data/uniport_pcp.rds")
pdb <- readRDS("~/PCP_model/data/pdb_pcp.rds")
all_pcp <- bind_rows(pdb,uniport)
all_pcp$pcp <- log2(all_pcp$pcp)
dup_seq <- all_pcp$seq[which(duplicated(all_pcp$seq))]
all_pcp <- all_pcp %>% filter(!(seq %in% dup_seq))
write.csv(all_pcp,file = "data/pcp_res_all.csv",row.names = F,quote = F)



