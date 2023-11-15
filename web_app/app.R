library(dplyr)
library(ggpubr)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(NGLVieweR)
library(waiter)
library(insight)
library(crayon)
library(bio3d)
library(shinyvalidate)
library(DT)

ifna <- function(x, elseval = NA) ifelse(is.na(x) || is.null(x), elseval, x)

cardb_seqs_adduser <- readRDS("/data/CAR-Shiny/cardb_seqs_adduser.rds")
###整理下数据库文件
cardb_seqs_adduser <- cardb_seqs_adduser %>% distinct_all()
user_file <- cardb_seqs_adduser %>% filter(grepl("^u",name))
non_user <- cardb_seqs_adduser %>% filter(!grepl("^u",name))
user_file <- user_file %>% filter(!(seq %in% non_user$seq))
cardb_seqs_adduser <- bind_rows(non_user,user_file)
saveRDS(cardb_seqs_adduser,file = "/data/CAR-Shiny/cardb_seqs_adduser.rds")

db_file <- readRDS("/data/CAR-Shiny/db_file.rds")

init_seq <- "DVVMTQTPLSLPVSLGDQASISCRSSQSLVHRNGNTYLHWYLQKPGQSPKLLIHKVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDLGVYFCSQSTHVPPLTFGAGTKLELKGGGGSGGGGSGGGGSEVQLLQSGPELEKPGASVMISCKASGSSFTGYNMNWVRQNIGKSLEWIGAIDPYYGGTSYNQKFKGRATLTVDKSSSTAYMHLKSLTSEDSAVYYCVSGMKYWGQGTSVTVSS"
cd19 <- "DIQMTQTTSSLSASLGDRVTISCRASQDISKYLNWYQQKPDGTVKLLIYHTSRLHSGVPSRFSGSGSGTDYSLTISNLEQEDIATYFCQQGNTLPYTFGGGTKLEITGSTSGSGKPGSGEGSTKGEVKLQESGPGLVAPSQSLSVTCTVSGVSLPDYGVSWIRQPPRKGLEWLGVIWGSETTYYNSALKSRLTIIKDNSKSQVFLKMNSLQTDDTAIYYCAKHYYYGGSYAMDYWGQGTSVTVSS"
header <- dashboardHeader(
  title = "CAR-PCP"
)

source("/data/CAR-Shiny/scripts/help_functions.R")

get_mut_pcp_pre <- function(pro_seq,car_fr,seq_id){

  mut_res <- get_mut_cars_all(car_fr = car_fr)
  mut_seq <- mut_res[[3]] %>% rename(seq=mut)
  write.table(mut_seq,
              file = paste0("/data/CAR-Shiny/tempdir/trans_pcp_",
                            seq_id,"_batch.txt"),
              row.names = F,quote = F,sep = "\t")
  
  car_fr$id <- 1:nrow(car_fr)
  intervals_k <- cumsum(c(1,nchar(car_fr$seq)))
  up_mut_pos <- as.numeric(mut_res[[1]])
  down_mut_pos <- as.numeric(mut_res[[2]])
  intervals_need_up <- findInterval(up_mut_pos,intervals_k)
  intervals_need_down <- findInterval(down_mut_pos,intervals_k)
  pos_need_up <- up_mut_pos - intervals_k[findInterval(up_mut_pos,intervals_k)] + 1
  pos_need_down <- down_mut_pos - intervals_k[findInterval(down_mut_pos,intervals_k)] + 1
  
  out_up <- car_fr %>% 
    rowwise() %>% 
    mutate(text_out=case_when(
      (grepl("FR",type) && (id %in% intervals_need_up)) ~ 
        ifelse(length(pos_need_up[intervals_need_up == id]) == 0,
               NA,
               highlight_aa(seq,pos_need_up[intervals_need_up == id],
                            type,"red")),
      (grepl("FR",type) && (!(id %in% intervals_need_up))) ~ 
        paste0("<span style='color:green;font-size:18px;'>",
               paste0(type,": ",seq),"</span>"),
      grepl("CDR",type) ~ paste0("<span style='color:blue;font-size:18px;'>",paste0(type,": ",seq),"</span>"),
      grepl("EXT",type) ~ paste0("<span style='color:black;font-size:18px;'>",paste0(type,": ",seq),"</span>")
    ))
  out_down <- car_fr %>% 
    rowwise() %>% 
    mutate(text_out=case_when(
      (grepl("FR",type) && (id %in% intervals_need_down)) ~ 
        ifelse(length(pos_need_down[intervals_need_down == id]) == 0,
               NA,
               highlight_aa(seq,pos_need_down[intervals_need_down == id],
                            type,"red")),
      (grepl("FR",type) && (!(id %in% intervals_need_down))) ~ 
        paste0("<span style='color:green;font-size:18px;'>",
               paste0(type,": ",seq),"</span>"),
      grepl("CDR",type) ~ paste0("<span style='color:blue;font-size:18px;'>",paste0(type,": ",seq),"</span>"),
      grepl("EXT",type) ~ paste0("<span style='color:black;font-size:18px;'>",paste0(type,": ",seq),"</span>")
    ))
  return(list(mut_seq,out_up,out_down))
}

# commd <- paste0("conda run -n dl ",
#                 "python /home/wt/CAR-Shiny/scripts/pre_pcp.py -i ",
#                 "/home/wt/CAR-Shiny/tempdir/trans_pcp_",seq_id,"_batch.txt ",
#                 "-o /home/wt/CAR-Shiny/tempdir/trans_pre_",seq_id,"_batch.txt ",
#                 "-m /home/wt/CAR-Shiny/pcp_tunning5/checkpoint-22176/")
# system(commd)

get_mut_plot <- function(pro_id,type){
  cardb_seqs_adduser <- readRDS("/data/CAR-Shiny/cardb_seqs_adduser.rds")
  ##type = q_k or k_q
  if (grepl("^u",pro_id)){
    mut_file <- readRDS(paste0("/data/CAR-Shiny/mut_cars/",pro_id,"_mut.rds"))
    mut_file$pcp <- NA
    for (i in 1:nrow(mut_file)){
      dt <- read.csv(paste0("/data/CAR-Shiny/mut_cars/user_res/",
                            mut_file$id[i],"_patches.csv"))
      mut_file$pcp[i] <- nrow(dt)
    }
    origin <- data.frame(pos="WT",
                         pcp=cardb_seqs_adduser$pcp[cardb_seqs_adduser$name==pro_id])
    mut <- bind_rows(mut_file %>% 
                       filter(grepl(substr(type,1,1),id)) %>% 
                       select(pos,pcp),
                     origin)
  }else{
    all_files <- readRDS(paste0("/data/CAR-Shiny/mut_cars/",type,"_file.rds")) %>% 
      mutate(origin_name = gsub(".rds","",origin_name))
    trans_name <- all_files$trans_name[all_files$origin_name==pro_id]
    files <- list.files(paste0("/data/CAR-Shiny/mut_cars/mut_seqs_",type,"_res/"),
                        pattern = 
                          paste0(trans_name, "_.+.csv"))
    res <- data.frame(id=gsub("_patches.csv","",files),
                      pcp = NA)
    for (i in 1:nrow(res)){
      dt <- read.csv(paste0("/data/CAR-Shiny/mut_cars/mut_seqs_",type,"_res/",
                            res$id[i],"_patches.csv"))
      res$pcp[i] <- nrow(dt)
    }
    mut <- readRDS(paste0("/data/CAR-Shiny/mut_cars/mut_seqs_",type,"/",
                          pro_id,".rds"))[[2]]
    mut$id <- paste0(trans_name,"_",c(1:nrow(mut)))
    mut <- left_join(mut,res)
    origin <- data.frame(pos="WT",
                         pcp=cardb_seqs_adduser$pcp[cardb_seqs_adduser$name==pro_id])
    mut <- bind_rows(mut %>% select(pos,pcp),
                     origin)
  }
  mut_plot <- mut %>% 
    mutate(mut_counts=lengths(strsplit(pos,split = ","))) %>% 
    mutate(mut_counts = ifelse(mut_counts==1,0,mut_counts)) %>% 
    arrange(mut_counts)
  mut_plot$pos <- factor(mut_plot$pos,levels = mut_plot$pos)
  mut_plot$`Mutation Counts` <- as.character(mut_plot$mut_counts)
  p <- ggbarplot(mut_plot, x = "pos", y = "pcp",
                 label = TRUE, label.pos = "out",
                 color = "Mutation Counts",fill="Mutation Counts",sort.val="desc")+
    rotate_x_text(45)+
    labs(x="Mutation Postion",y="PCP score")
  res <- list(mut=mut,p=p)
  return(res)
}

get_mut_plot_model <- function(pre_res,type){
  mut <- pre_res %>% 
    filter(grepl(type,id))
  mut_plot <- mut %>% 
    mutate(mut_counts=lengths(strsplit(pos,split = ","))) %>% 
    mutate(mut_counts = ifelse(mut_counts==1,0,mut_counts)) %>% 
    arrange(mut_counts)
  mut_plot$pos <- factor(mut_plot$pos,levels = mut_plot$pos)
  mut_plot$`Mutation Counts` <- as.character(mut_plot$mut_counts)
  p <- ggbarplot(mut_plot, x = "pos", y = "pre",
                 label = TRUE, label.pos = "out",
                 color = "Mutation Counts",fill="Mutation Counts",sort.val="desc")+
    rotate_x_text(45)+
    labs(x="Mutation Postion",y="PCP score")
  res <- list(mut=mut,p=p)
}

get_mut_highlight <- function(pro_id){
  cdr_fr <- read.csv(paste0("/data/CAR-Shiny/cdr_fr/",pro_id,".csv"))
  cdr_fr$id <- 1:nrow(cdr_fr)
  mut_car_k <- readRDS(paste0("/data/CAR-Shiny/mut_cars/mut_seqs_k_q/",pro_id,".rds"))
  intervals_k <- cumsum(c(1,nchar(cdr_fr$seq)))
  intervals_need_k <- findInterval(mut_car_k$top5_mut,intervals_k)
  pos_need_k <- mut_car_k$top5_mut -intervals_k[findInterval(mut_car_k$top5_mut,intervals_k)] + 1
  
  mut_car_q <- readRDS(paste0("/data/CAR-Shiny/mut_cars/mut_seqs_q_k/",pro_id,".rds"))
  intervals_q <- cumsum(c(1,nchar(cdr_fr$seq)))
  intervals_need_q <- findInterval(mut_car_q$top5_mut,intervals_q)
  pos_need_q <- mut_car_q$top5_mut -intervals_q[findInterval(mut_car_q$top5_mut,intervals_q)] + 1
  
  intervals_need <- c(intervals_need_k,intervals_need_q)
  pos_need <- c(pos_need_k,pos_need_q)
  cdr_fr <- cdr_fr %>% 
    rowwise() %>% 
    mutate(text_out=case_when(
      (grepl("FR",type) && (id %in% intervals_need)) ~ 
        ifelse(length(pos_need[intervals_need == id]) == 0,
               NA,
               highlight_aa(seq,pos_need[intervals_need == id],
                            type)),
      (grepl("FR",type) && (!(id %in% intervals_need))) ~ 
        paste0("<span style='color:green;font-size:18px;'>",
               paste0(type,": ",seq),"</span>"),
      grepl("CDR",type) ~ paste0("<span style='color:blue;font-size:18px;'>",paste0(type,": ",seq),"</span>"),
      grepl("EXT",type) ~ paste0("<span style='color:black;font-size:18px;'>",paste0(type,": ",seq),"</span>")
    ))
  return(cdr_fr)
}

get_fa_from_seq <- function(seq){
  seqs <- strsplit(seq,">")[[1]]
  seqs <- seqs[2:length(seqs)]
  res <- vector("list",length = length(seqs))
  for (i in seq_along(res)){
    name_seq <- strsplit(seqs[i],"\n")[[1]]
    names(res)[i] <- name_seq[1]
    res[[i]] <- name_seq[2]
  }
  return(res)
}

fa2df <- function(fa){
  seqs <- sapply(fa,as.character)
  names <- names(fa)
  dt <- data.frame(
    ids = names,
    seq = seqs
  )
  return(dt)
}


sidebar <- dashboardSidebar(
  sidebarMenu(id = "tabs",
    #menuItem("CARDb", tabName = "cardb"),
    # menuItem("PCP-Calculator", tabName = "pcp",
    #          menuSubItem("Batch calculator","batch_tab"),
    #          menuSubItem("Tunning","tunning_tab"))
    menuItem("Single calculator", tabName = "tunning_tab"),
    menuItem("Batch calculator", tabName = "batch_tab")
  )
)

body <- dashboardBody(
  useWaiter(),
  tabItems(
    # tabItem(tabName = "cardb",
    #         fluidRow(
    #           box(
    #             title = "Input", status = "primary", solidHeader = TRUE,
    #             collapsible = FALSE,width=4,
    #             ##选择库里面的 CAR
    #             pickerInput(
    #               inputId = "cardb_id",
    #               label = "Select one CAR", 
    #               choices = db_file$name,
    #               options = list(
    #                 style = "material-flat")
    #             )
    #           ),
    #           valueBoxOutput(outputId = "cardb_pcp_score",width=3)
    #         ),
    #         fluidRow(
    #           box(
    #             title = "PDB View", status = "primary", solidHeader = TRUE,
    #             collapsible = FALSE,width=4,
    #             ##展示 PDB 
    #             shinycssloaders::withSpinner(
    #               NGLVieweROutput("cardb_pdb") 
    #             )
    #           ),
    #           box(
    #             title = "scFv Sequence", status = "primary", solidHeader = TRUE,
    #             collapsible = FALSE,width=6,
    #             ###展示 CDR,FR 结构，并展示前 5 个突变位点
    #             shinycssloaders::withSpinner(
    #               htmlOutput("cardb_seq")
    #             )
    #           )
    #         ),
    #         fluidRow(
    #           box(
    #             title = "Tuning up PCP (Q→K)", status = "primary", solidHeader = TRUE,
    #             collapsible = FALSE,width=5,
    #             ##展示柱状图，Q-K 不同突变的 PCP 值
    #             shinycssloaders::withSpinner(
    #               plotOutput("cardb_up")
    #             )
    #           ),
    #           box(
    #             title = "Tuning down PCP (K→Q)", status = "primary", solidHeader = TRUE,
    #             collapsible = FALSE,width=5,
    #             ##展示柱状图，K-Q 不同突变的 PCP 值
    #             shinycssloaders::withSpinner(
    #               plotOutput("cardb_down")
    #             )
    #           )
    #         )
    # ),
    
    tabItem(tabName = "batch_tab",
            bootstrapPage(
              ###input
              fluidRow(
                box(title = "Input",width=40, solidHeader = FALSE, collapsible = FALSE,
                    fluidRow(
                      ##输入序列
                      column(width=8,
                             textAreaInput("batch_seq",
                                           label = "Input scFV sequences (In FASTA format)", 
                                           value = paste0(">GD2\n",init_seq,"\n",
                                                          ">CD19\n",cd19),
                                           height="140px"
                                           )
                      )
                    ),
                    fluidRow(
                      ##输入文件
                      column(width=8,
                             fileInput("fasta_input", "Choose FASTA file",
                                       multiple = FALSE,
                                       accept = ".fasta")
                      )
                    ),
                    fluidRow(
                      column(width = 12,
                             actionBttn(
                               inputId = "run_batch",
                               label = "Run",
                               color = "primary",
                               style = "gradient",
                               block = FALSE,size="md"
                             )
                      )
                    )
                )
              ),
              ##output
              fluidRow(
                uiOutput("batch_out")##输出一个表格
              )
            )
      
    ),
    
    tabItem(tabName = "tunning_tab",
            bootstrapPage(
              
              fluidRow(
                box(title = "Input",width=8, solidHeader = FALSE, 
                    collapsible = FALSE,
                    fluidRow(
                      column(width = 12,
                             textInput("seq",
                                       label = "Input sequence", 
                                       value = init_seq)
                      )
                    ),
                    fluidRow(
                      column(width = 3,
                             actionBttn(
                               inputId = "run_tunning",
                               label = "Run",
                               color = "primary",
                               style = "gradient",
                               block = TRUE,size="md"
                             )),
                      column(width = 3,
                             uiOutput("ready_tunning"))
                    )
                ),
                valueBoxOutput("s_pcp_score", width = 2)
              ),
              fluidRow(
                valueBoxOutput("mut_seq_num", width = 3),
                valueBoxOutput("es_time", width = 3)
              ),
              fluidRow(
                uiOutput("opt_seq_up_ui"),
                uiOutput("opt_seq_down_ui")
              ),
              fluidRow(
                uiOutput("tunning_res_ui"),
                uiOutput("up_or_down_res_ui")
              )
            )
    )
  )
)

ui <- dashboardPage(header, sidebar, body)

server <- function(input, output) {
  w <- Waiter$new(
    html = spin_3k(), 
    color = transparent(.5)
  )
  ###数据库模块
  # output$cardb_pcp_score <- renderValueBox({
  #   db_file <- readRDS("/data/CAR-Shiny/db_file.rds")
  #   valueBox(
  #     "PCP-Score", 
  #     tags$p(db_file$pcp[which(db_file$name==input$cardb_id)],
  #            style = "font-size: 150%;"), 
  #     icon = icon("calculator"),
  #     color = "yellow"
  #   )
  # })
  # output$cardb_pdb <- renderNGLVieweR({
  #   patches <- read.csv(paste0("/data/CAR-Shiny/local_pipeline_res/",
  #                              input$cardb_id,"_patches.csv"))
  #   patch1 <- paste(patches$residue_number[which(patches$patch=="p1")],
  #                   collapse = " or ")
  #   patch2 <- paste(patches$residue_number[which(patches$patch=="p2")],
  #                   collapse = " or ")
  #   patch3 <- paste(patches$residue_number[which(patches$patch=="p3")],
  #                   collapse = " or ")
  #   NGLVieweR(paste0("/data/CAR-Shiny/local_pipeline_res/",input$cardb_id,"_model.pdb")) %>% 
  #     addRepresentation("surface",
  #                       param = list(
  #                         colorValue = "white"
  #                       )
  #     ) %>%  addRepresentation("surface",
  #                              param = list(
  #                                sele = patch1,
  #                                labelType = "format",
  #                                labelFormat = "[%(resname)s]%(resno)s", # or enter custom text
  #                                labelGrouping = "residues", # or "atom" (eg. sele = "20:A.CB")
  #                                color = "#0000C6",
  #                                fontFamiliy = "sans-serif",
  #                                xOffset = 1,
  #                                yOffset = 0,
  #                                zOffset = 0,
  #                                fixedSize = TRUE,
  #                                radiusType = 1,
  #                                radiusSize = 2, # Label size
  #                                showBackground = FALSE
  #                                # backgroundColor="black",
  #                                # backgroundOpacity=0.5
  #                              )
  #     ) %>% 
  #     addRepresentation("surface",
  #                       param = list(
  #                         sele = patch2,
  #                         labelType = "format",
  #                         labelFormat = "[%(resname)s]%(resno)s", # or enter custom text
  #                         labelGrouping = "residues", # or "atom" (eg. sele = "20:A.CB")
  #                         color = "#1383FD",
  #                         fontFamiliy = "sans-serif",
  #                         xOffset = 1,
  #                         yOffset = 0,
  #                         zOffset = 0,
  #                         fixedSize = TRUE,
  #                         radiusType = 1,
  #                         radiusSize = 2, # Label size
  #                         showBackground = FALSE
  #                         # backgroundColor="black",
  #                         # backgroundOpacity=0.5
  #                       )
  #     ) %>% 
  #     addRepresentation("surface",
  #                       param = list(
  #                         sele = patch3,
  #                         labelType = "format",
  #                         labelFormat = "[%(resname)s]%(resno)s", # or enter custom text
  #                         labelGrouping = "residues", # or "atom" (eg. sele = "20:A.CB")
  #                         color = "#9CCBFE",
  #                         fontFamiliy = "sans-serif",
  #                         xOffset = 1,
  #                         yOffset = 0,
  #                         zOffset = 0,
  #                         fixedSize = TRUE,
  #                         radiusType = 1,
  #                         radiusSize = 2, # Label size
  #                         showBackground = FALSE
  #                         # backgroundColor="black",
  #                         # backgroundOpacity=0.5
  #                       )
  #     )
  #   
  # })
  # 
  # output$cardb_seq <- renderUI({
  #   cdr_fr <- get_mut_highlight(pro_id = input$cardb_id)
  #   HTML(paste(cdr_fr$text_out, collapse = '<br/>'))
  # })
  # 
  # output$cardb_up <- renderPlot({
  #   res <- get_mut_plot(pro_id = input$cardb_id,type = "q_k")
  #   res$p
  # })
  # 
  # output$cardb_down <- renderPlot({
  #   res <- get_mut_plot(pro_id = input$cardb_id,type = "k_q")
  #   res$p
  # })
  
  ###计算 PCP
  ###batch mode
  batch_seq_pcp <- eventReactive(input$run_batch,{
    seq_id <- paste0("u", gsub("[.].+","",as.character(as.numeric(Sys.time()))))
    if (!is.null(input$fasta_input)){
      dt <- seqinr::read.fasta(input$fasta_input$datapath,
                               seqtype = "AA",as.string = TRUE)
    }else{
      seqs <- get_fa_from_seq(input$batch_seq)
      seqinr::write.fasta(sequences = seqs,
                          names = names(seqs),
                          file.out = paste0("/data/CAR-Shiny/tempdir/",seq_id,"_batch.fa"))
      dt <- seqinr::read.fasta(paste0("/data/CAR-Shiny/tempdir/",seq_id,"_batch.fa"),
                               seqtype = "AA",as.string = TRUE)
    }
    ##计算
    dt_seqs <- fa2df(dt)
    write.table(dt_seqs,
                file = paste0("/data/CAR-Shiny/tempdir/trans_pcp_",
                                 seq_id,"_batch.txt"),
                row.names = F,quote = F,sep = "\t")
    commd <- paste0("conda run -n dl ",
                    "python /data/CAR-Shiny/scripts/pre_pcp.py -i ",
                    "/data/CAR-Shiny/tempdir/trans_pcp_",seq_id,"_batch.txt ",
                    "-o /data/CAR-Shiny/tempdir/trans_pre_",seq_id,"_batch.txt ",
                    "-m /data/CAR-Shiny/pcp_tunning5/checkpoint-22176/")
    system(commd)
    dt_pre <- read.table(paste0("/data/CAR-Shiny/tempdir/trans_pre_",
                                seq_id,"_batch.txt"))
    dt_pre$pre <- round(2^(dt_pre$pre))
    return(dt_pre)
  })
  
  observeEvent(input$run_batch,{
    output$batch_out <- renderUI({
      box(title = "PCP Score", status = "primary", solidHeader = TRUE,width=10,
          shinycssloaders::withSpinner(DTOutput("batch_res"))
      )
    }) 
  })
  
  output$batch_res <- DT::renderDataTable(
    DT::datatable(batch_seq_pcp() %>% 
                    rename(ID = ids,
                           Sequence = seq,
                           `PCP Score` = pre),
		  rownames = as.character(1:nrow(batch_seq_pcp()))) %>% 
      DT::formatStyle(
        columns = "Sequence", 
        display = "block", 
        `overflow-wrap`= 'break-word',
        overflow = 'hidden',
        width = "780px"
      ),
    options = list(lengthChange = TRUE,pageLength = 10), selection = 'none'
  )
  ###tunning mode
  ##UI
  pcp_pre <- eventReactive(input$run_tunning,{
    cardb_seqs_adduser <- readRDS("/data/CAR-Shiny/cardb_seqs_adduser.rds")
    if (input$seq %in% cardb_seqs_adduser$seq){
      seq_id <- cardb_seqs_adduser$name[which(cardb_seqs_adduser$seq == input$seq)]
      seq_id <- seq_id[1]
    }else{
      seq_id <- paste0("u", gsub("[.].+","",as.character(as.numeric(Sys.time()))))
    }
    ##输入序列，得到dataframe
    temp_file <- paste0("/data/CAR-Shiny/tempdir/",seq_id,".fa")
    seqinr::write.fasta(input$seq,names = seq_id,file.out = temp_file)
    commd <- paste0("/data/CAR-Shiny/AbRSA/AbRSA -i ", temp_file, 
                    " -g -o /data/CAR-Shiny/tempdir/temp_number_",seq_id,
                    " > /data/CAR-Shiny/tempdir/","out_",seq_id)
    system(commd)
    out <- data.table::fread(paste0("/data/CAR-Shiny/tempdir/out_",seq_id),
                             skip = "#similarity",
                             data.table = F,sep = ":",header = F)
    colnames(out) <- c("type","seq")
    res <- get_mut_pcp_pre(input$seq,out,seq_id)
    
    ###自身的 PCP
    dt <- data.frame(seq=input$seq)
    write.table(dt,
                file = paste0("/data/CAR-Shiny/tempdir/trans_pcp_",
                              seq_id,"_wt.txt"),
                row.names = F,quote = F,sep = "\t")
    commd <- paste0("source ~/.bashrc;conda run -n dl ",
                    "python /data/CAR-Shiny/scripts/pre_pcp.py -i ",
                    "/data/CAR-Shiny/tempdir/trans_pcp_",seq_id,"_wt.txt ",
                    "-o /data/CAR-Shiny/tempdir/trans_pre_",seq_id,"_wt.txt ",
                    "-m /data/CAR-Shiny/pcp_tunning5/checkpoint-22176/")
    system(commd)
    dt_pre <- read.table(paste0("/data/CAR-Shiny/tempdir/trans_pre_",
                                seq_id,"_wt.txt"))
    dt_pre$pre <- round(2^(dt_pre$pre))
    res[[4]] <- dt_pre$pre
    if (!(input$seq %in% cardb_seqs_adduser$seq)){
      new_seq <- data.frame(name=seq_id,seq=input$seq,pcp=res[[4]])
      cardb_seqs_adduser <- bind_rows(cardb_seqs_adduser,new_seq)
      saveRDS(cardb_seqs_adduser,file = "/data/CAR-Shiny/cardb_seqs_adduser.rds")
    }
    return(res)
  })
  
  pcp_cal <- eventReactive(input$run_tunning_check,{
    res <- pcp_pre()
    cardb_seqs_adduser <- readRDS("/data/CAR-Shiny/cardb_seqs_adduser.rds")
    seq_id <- cardb_seqs_adduser$name[which(cardb_seqs_adduser$seq == input$seq)]
    seq_id <- seq_id[1]
    if (!file.exists(paste0("/data/CAR-Shiny/tempdir/trans_pre_",
                            seq_id,"_batch.txt"))){
      commd <- paste0("source ~/.bashrc;conda run -n dl ",
                      "python /data/CAR-Shiny/scripts/pre_pcp.py -i ",
                      "/data/CAR-Shiny/tempdir/trans_pcp_",seq_id,"_batch.txt ",
                      "-o /data/CAR-Shiny/tempdir/trans_pre_",seq_id,"_batch.txt ",
                      "-m /data/CAR-Shiny/pcp_tunning5/checkpoint-22176/")
      system(commd)
    }
    dt_pre <- read.table(paste0("/data/CAR-Shiny/tempdir/trans_pre_",
                                seq_id,"_batch.txt"))
    dt_pre$pre <- round(2^(dt_pre$pre))
    return(dt_pre)
  })
  
  # get_plot <- eventReactive(input$run_tunning_check,{
  #   pcp_pre_res <- pcp_pre()
  #   pcp_cal_res <- pcp_cal()
  #   if (pcp_pre_res[[4]] > 56){
  #     res <- pcp_cal_res %>% 
  #       filter(pre>=46 & pre<=56) %>% 
  #       arrange(pre) %>% 
  #       slice_head(n=15)
  #     if (nrow(res) == 0){
  #       res <- pcp_cal_res %>% 
  #         mutate(diff = pre-46) %>% 
  #         filter(pre <= 46) %>% 
  #         arrange(desc(diff)) %>% 
  #         slice_head(n=15)
  #     }
  #   } else if (pcp_pre_res[[4]] < 46){
  #     res <- pcp_cal_res %>% 
  #       filter(pre>=46 & pre<=56) %>% 
  #       arrange(pre) %>% 
  #       slice_head(n=15)
  #     if (nrow(res) == 0){
  #       res <- pcp_cal_res %>% 
  #         filter(pre > pcp_pre_res[[4]]) %>% 
  #         arrange(pre) %>% 
  #         slice_head(n=15)
  #     }
  #   }else{
  #     res <- pcp_cal_res %>% 
  #       filter(pre>=46 & pre<=56) %>% 
  #       arrange(pre) %>% 
  #       slice_head(n=15)
  #   }
  #   res <- res %>% 
  #     select(pos,origin_aa,mut_aa,pre) %>% 
  #     mutate(x = paste(pos,origin_aa,mut_aa,sep  = "-"))
  #   res <- res %>% 
  #     mutate(mut_counts=lengths(strsplit(pos,split = ","))) %>% 
  #     mutate(mut_counts = ifelse(mut_counts==1,0,mut_counts))
  #   res$`Mutation Counts` <- as.character(res$mut_counts)
  #   p <- ggbarplot(res, x = "x", y = "pre",
  #                  label = TRUE, label.pos = "out",
  #                  color = "Mutation Counts",fill="Mutation Counts",sort.val="desc")+
  #     rotate_x_text(45)+
  #     labs(x="Mutation Postion",y="PCP score")
  #   return(p)
  # })
  
  observeEvent(input$run_tunning,{
    w$show()
    res <- pcp_pre()
    w$hide()
    # output$overview_ui <- renderUI({
    #   box(title = "Overview", status = "primary", 
    #       solidHeader = TRUE,
    #       collapsible = FALSE,width=4,
    #       shinycssloaders::withSpinner(
    #         htmlOutput("overview") 
    #       ),
    #       actionBttn(
    #         inputId = "run_tunning_check",
    #         label = "Tunning",
    #         color = "primary",
    #         style = "gradient",
    #         block = FALSE,size="md"
    #       )
    #   )
    # })
    output$s_pcp_ui <- renderUI({
      shinycssloaders::withSpinner(valueBoxOutput("s_pcp_score", width = 3))
    })
    output$s_pcp_score <- renderValueBox({
      res <- pcp_pre()
      valueBox(
        res[[4]],"PCP Score",icon = icon("calculator",lib = "font-awesome"),
        color = "yellow",width = 3
      )
    })
    output$mut_seq_num <- renderValueBox({
      res <- pcp_pre()
      valueBox(
        nrow(res[[1]]), "Number of Mutated Sequences (1-4)", 
        icon = icon("hashtag",lib = "font-awesome"),
        color = "purple",width = 3
      )
    })
    output$es_time <- renderValueBox({
      res <- pcp_pre()
      valueBox(
        round(nrow(res[[1]])/102), "Estimated Time (Minutes)", icon = icon("clock",lib = "font-awesome"),
        color = "purple",width = 3
      )
    })
    output$ready_tunning <- renderUI({
      actionBttn(
        inputId = "run_tunning_check",
        label = "Tunning",
        color = "primary",
        style = "gradient",
        block = TRUE,size="md"
      )
    })
    
    output$opt_seq_up_ui <- renderUI({
      box(
        title = "Tunning Up (Q → K)", status = "primary", solidHeader = TRUE,
        collapsible = FALSE,width=6,
        shinycssloaders::withSpinner(
          htmlOutput("opt_seq_up")
        )
      )
    })
    output$opt_seq_down_ui <- renderUI({
      box(
        title = "Tunning Down (K → Q)", status = "primary", solidHeader = TRUE,
        collapsible = FALSE,width=6,
        shinycssloaders::withSpinner(
          htmlOutput("opt_seq_down")
        )
      )
    })
  })
  
  observeEvent(input$run_tunning_check,{
    w$show()
    res <- pcp_cal()
    w$hide()
    output$up_or_down_res_ui <- renderUI({
      box(
        title = "PCP Score Distribution", 
        status = "primary", solidHeader = TRUE,
        collapsible = FALSE,width=10,
        shinycssloaders::withSpinner(
          plotOutput("up_or_down_res",
                     height = 500, 
                     brush = brushOpts(id = "plot_brush", direction = "x",
                                       fill = "blue", opacity = 0.5))
        ),
        DTOutput("filter_res")
      )
    })
    output$tunning_res_ui <- renderUI({
      box(
        title = "All Tunning Results", 
        status = "primary", solidHeader = TRUE,
        collapsible = FALSE,width=10,
        fluidRow(
          column(width = 10,
                 shinycssloaders::withSpinner(
                   DTOutput("tunning_res")
                   )
                 ),
          column(width = 2,
                 shinycssloaders::withSpinner(
                   downloadLink('download_all_tunning',"Download")
                   )
                 )
        )
      )
    })
  })
  ##cal
  output$overview <- renderUI({
    res <- pcp_pre()
    show_text <- paste0('<p style="font-size:16px; "><b>PCP score:</b> ',
                        res[[4]],"</p>",
                        '<p style="font-size:16px; "><b>Number of mutated sequences:</b> ',
                        nrow(res[[1]]),"</p>",
                        '<p style="font-size:16px; "><b>Estimated time:</b> ',
                        round(nrow(res[[1]])/102),
                        " min</p>")
    HTML(show_text)
  })
  output$opt_seq_up <- renderUI({
    res <- pcp_pre()
    HTML(paste(res[[2]]$text_out,collapse = '<br/>'))
  })
  output$opt_seq_down <- renderUI({
    res <- pcp_pre()
    HTML(paste(res[[3]]$text_out,collapse = '<br/>'))
  })
  
  output$tunning_res <- DT::renderDataTable(
    DT::datatable(pcp_cal() %>% 
                    select(-c("seq","origin_seq","ll")) %>% 
                    rename(`Mutation Position`=pos,
                           `WT AA`=origin_aa,
                           `Mut AA`=mut_aa,
                           `Tunning Type`=tunning_type,
                           `PCP Score`=pre)),
    options = list(lengthChange = TRUE,pageLength = 8), 
    selection = 'none'
  )
  output$download_all_tunning <- downloadHandler(
      filename = function() {
        paste('data-', Sys.Date(), '.xlsx', sep='')
      },
      content = function(con) {
        dt <- pcp_cal() %>% select(-ll) %>% 
          rename(mut_seq = seq)
        WriteXLS::WriteXLS("dt",ExcelFileName = con)
      }
    )
  output$up_or_down_res <- renderPlot({
    pcp_cal_res <- pcp_cal()
    ggdensity(data=pcp_cal_res,x="pre",fill = "#00AFBB",xlab = "PCP Score",
              ylab = "Density")+
      geom_vline(xintercept = 46, 
                 color = "red", linewidth=2)+
      geom_vline(xintercept = 56,
                 color = "red", linewidth=2)
  })
  selmin <- reactive(ifna(input$plot_brush$xmin, elseval = 46))
  selmax <- reactive(ifna(input$plot_brush$xmax, elseval = 56))
  output$filter_res <- DT::renderDataTable(
    DT::datatable(pcp_cal() %>% 
                    select(-c("seq","origin_seq","ll")) %>% 
                    filter(pre >= selmin() & pre <= selmax()) %>% 
                    rename(`Mutation Position`=pos,
                           `WT AA`=origin_aa,
                           `Mut AA`=mut_aa,
                           `Tunning Type`=tunning_type,
                           `PCP Score`=pre)),
    options = list(lengthChange = TRUE,pageLength = 8), 
    selection = 'none'
  )
}

shinyApp(ui, server)
