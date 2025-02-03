library("RamanEx") #
func_lib()
library('uwot')
library("plotly")
library("paletteer")

img_format <- c("svg","pdf")

folderpath <- "---"
outpath <- func_initialization(folderpath, output_folder = "output")
Name_group <- c("group_strain", "conc", "time") #, "group_conc", "group_time"

read_allsinspc(folderpath, Name_group,outpath,instrument = "X")

##################### 2. read and reorganize
filepath <- paste(folderpath,"/output/files",sep = "")
all_data_list <- read_3data(filepath)

spc_reor_filt <- dcast(all_data_list[[1]],filename ~ wavenumber)
meta_filt <- filter(all_data_list[[2]],filename %in% spc_reor_filt$filename) #,wavenum == "854", group_strain %in%c("JP1","JP2")
spc_reor_filt <- spc_reor_filt %>% filter(filename %in% meta_filt$filename)
data_info <- all_data_list[[3]]

spc_reor_filt <- pre_reorganize(spc_reor_filt,meta_filt)

##################### 3. data clean & QC
spc_reor_filt_qc <- cal_QC(spc_reor_filt)

spc_qc_good <- spc_reor_filt_qc %>% filter(Pre_Noise_mean < 0.1, Pre_Noise_sd < 0.1,
                                           Pre_Signal_max > 0.1, Pre_Signal_max < 1.2,
                                           Ori_Intensity_max < 65500,Ori_Intensity_min >= -2500) #

spc_good <- spc_reor_filt %>% filter(filename %in% spc_qc_good$filename)
hs_ori <- unique(spc_good[,c(TRUE,!duplicated(as.numeric(colnames(spc_reor_filt[,-1]))))])

hs_meta <- meta_filt %>% filter(filename %in% hs_ori$filename)
hs_ori <- hs_ori %>% filter(filename %in% hs_meta$filename)

##################### 4. pretreatment
hs_baseline <- pre_baseline(hs_ori)
hs_norm <- pre_norm(hs_baseline,method = "max")
hs_pre <- pre_spike_matrix(hs_norm)

outpath_files <- func_initialization(outpath, output_folder = "files")
fwrite(hs_pre, paste(outpath_files, "/alldata_pre", ".txt", sep = ""),
       row.names = F, col.names = T, quote = F, sep = "\t")


hs_CD_ratio <- cal_CDR(hs_pre,hs_meta)
fwrite(hs_CD_ratio, paste(outpath_files, "/CDR_groups", ".txt", sep = ""),
       row.names = F, col.names = T, quote = F, sep = "\t")
# ################### 5. UMAP
outpath_files <- func_initialization(outpath, output_folder = "files")

hs_CD_ratio <- fread(paste(outpath_files, "/CDR_groups.txt", sep = ""), header = TRUE, sep = "\t")
hs_pre <- fread(paste(outpath_files, "/alldata_pre.txt", sep = ""), header = TRUE, sep = "\t")
hs_meta_1 <- hs_CD_ratio %>% filter(group_strain %in% c("S1","S2","S3","S4"))
hs_pre_1 <- hs_pre %>% filter(filename %in% hs_meta_1$filename)

set.seed(50) #44 

data.umap <- uwot::umap(hs_pre_1[,-1])#spc_waveumber(hs_pre_1,400,1750)[,-1]
data.umap <- as.data.frame(data.umap)
colnames(data.umap) <- c('UMAP_1','UMAP_2')
#data.umap$UMAP_2 <- data.umap$UMAP_2 * -1
data.umap <- cbind(data.umap,filename = hs_pre_1$filename)
data.umap <- left_join(data.umap,hs_meta_1, by = "filename")
data.umap <- data.umap %>% filter(CD_ratio > -0.05)

fwrite(data.umap, paste(outpath_files, "/data-umap", ".txt", sep = ""),
       row.names = F, col.names = T, quote = F, sep = "\t")

data.umap <- fread(paste(outpath_files, "/data-umap.txt", sep = ""), header = TRUE, sep = "\t")

color_hetero_props <- c("#1f4047","#1e592b","#8c9976","#adbd95","#c4d18c",
                        "#f1cd9b","#f6b654","#f0634a","#dc3e32","#a73436")
color_YLGn_props <- colorRampPalette(rev(brewer.pal(1,"YlGn")))
data.umap$CD_ratio <- as.numeric(as.character(data.umap$CD_ratio))
inds <- data.umap_FQ$filename

data.umap_FQ <- data.umap %>% filter(group_strain %in% c("S1","S2","S3","S4"))
data.umap_FQ$filename <- factor(data.umap_FQ$filename,levels = inds)
umap_plot <- ggplot(data.umap_FQ, aes(UMAP_1, UMAP_2,group = group_strain,color = CD_ratio)) + 
  geom_point(size=1) + theme_bw() + 
  scale_color_gradientn(na.value = "#D8E5E9",colors = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_all-0424.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)


# p<-plot_ly(data.umap_FQ,
#            x=~UMAP_1,
#            y=~UMAP_2,
#            text=~row.names(data.umap_FQ),
#            color=~CD_ratio)
# p
# htmlwidgets::saveWidget(as_widget(p), "graphU.html")



#### FigA_2 
color_hetero_props <- c("#1f4047","#1e592b","#8c9976","#adbd95","#c4d18c",
                        "#f1cd9b","#f6b654","#f0634a","#dc3e32","#a73436")

data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S2"|
                                data.umap_FQ_1$group_strain == "S3"|
                                data.umap_FQ_1$group_strain == "S4")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_1-0424.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)


#### FigA_2 
data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S1"|
                                data.umap_FQ_1$group_strain == "S3"|
                                data.umap_FQ_1$group_strain == "S4")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_2-0424.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)

data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S1"|
                                data.umap_FQ_1$group_strain == "S2"|
                                data.umap_FQ_1$group_strain == "S4")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_3-0424.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)

data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S1"|
                                data.umap_FQ_1$group_strain == "S2"|
                                data.umap_FQ_1$group_strain == "S3")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_4-0424.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)

data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S2"|data.umap_FQ_1$group_strain == "S3")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_14.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)

data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S2"|data.umap_FQ_1$group_strain == "S4")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_13.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)

data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S3"|data.umap_FQ_1$group_strain == "S4")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_12.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)

data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S1"|data.umap_FQ_1$group_strain == "S2")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_34.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)


data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S1"|data.umap_FQ_1$group_strain == "S4")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_23.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)


data.umap_FQ_1 <- data.umap_FQ
data.umap_FQ_1$CD_ratio[which(data.umap_FQ_1$group_strain == "S1"|data.umap_FQ_1$group_strain == "S3")] <- NA
data.umap_FQ_1 <- data.umap_FQ_1[sort(order(data.umap_FQ_1$CD_ratio,decreasing = TRUE),decreasing = TRUE),]
data.umap_FQ_1$filename <- factor(data.umap_FQ_1$filename,levels = inds)

umap_plot <- ggplot(data.umap_FQ_1, aes(UMAP_1, UMAP_2,color = CD_ratio)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.8) + theme_bw() + 
  scale_colour_gradientn(na.value = "#D8E5E9",colours = color_hetero_props,guide = "colourbar",limits = c(0,max(data.umap_FQ$CD_ratio))) +
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_FQ_24.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 7,height = 6)

data.umap_JN <- data.umap %>% filter(group_strain %in% c("JP1","JP2"))

umap_plot <- ggplot(data.umap_JN, aes(UMAP_1, UMAP_2,color = group_strain)) + 
  geom_point(aes(color = CD_ratio),size=1.5,alpha = 0.5) + theme_bw() + 
  scale_colour_gradientn(colours = color_hetero_props,guide = "colourbar") +
  facet_grid(.~group_strain) + 
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_JN.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 10,height = 6)


umap_plot <- ggplot(data.umap_JN, aes(UMAP_1, UMAP_2,color = group_strain)) + 
  geom_point(aes(color = group_strain),size=1.5,alpha = 0.3) + theme_bw() +
  scale_fill_manual(values = c("#336699","#ED550A"),aesthetics = c("fill","color")) +
  facet_grid(.~wavenum) + 
  theme(panel.grid = element_blank(), 
        panel.grid.minor = element_blank(), 
        # legend.title = element_blank(), 
        # legend.text = element_markdown(size = 15), 
        # legend.background = element_blank(), text = element_text(color = "black"), 
        axis.text.x = element_text(size = 15, angle = 0), 
        axis.text.y = element_text(size = 15), 
        axis.ticks = element_line(linewidth  = 1), 
        axis.ticks.length = unit(0.4, "lines"), 
        axis.title = element_text(size = 15))
umap_plot
ggsave(filename=paste(outpath,"/","plotpoint_JN-3.svg", sep=""),plot=umap_plot,
       limitsize=F,width = 6,height = 6)



#### 特征峰
hs_pre_melt <- melt(hs_pre_1,id.vars = c("filename"),variable.name = "wavenumber",value.name = "value")
RICH_data_wavenumber <- as.numeric(as.character(sort(unique(hs_pre_melt$wavenumber))))

portein <- c(1004,1619)
nucleic <- c(786,814) 
lipid <- c(719, 1082)
alcohol <- c(877,1046)

wave_all <- portein

RICH_peaks <- c()
for( i in wave_all){
  
  RICH_peaks <- c(RICH_peaks,RICH_data_wavenumber[which.min(abs(RICH_data_wavenumber - i))])
}

RICH_slect_spc <- filter(hs_pre_melt,wavenumber %in% RICH_peaks)
RICH_slect_spc <- left_join(RICH_slect_spc,hs_CD_ratio,by = "filename")

plot <-ggplot(RICH_slect_spc,aes(x = factor(group_strain),y = value,group = group_strain,fill = group_strain)) +
  geom_point(aes(group = group_strain, color = group_strain), position = position_jitter(0.15),alpha = 0.2) +
  geom_boxplot(aes(group = group_strain),alpha = 0.2,outlier.alpha = 0) +
  scale_fill_manual(values = c("#FBA240","#104D68","#FA8300",
                               "#1F799F","#BF7930","#61B4CF"),aesthetics = c("fill","color")) +
  facet_grid(. ~ wavenumber) +
  default_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 15, angle = 45,hjust = 1,vjust = 1)) +
  xlab("Group") +
  ylab('Metabolic activity (CD-ratio)')
plot
ggsave(filename=paste(outpath,"/","plot_box_allwave_","portein",".svg", sep=""),plot=plot,
       limitsize=T,width = 10,height = 6)

hist_CDR_plot <- ggplot(RICH_slect_spc, aes(x = value, fill = factor(group_strain))) +
  geom_density(aes(y = ..ndensity..),position = "identity",alpha = 0.4) + 
  geom_histogram(aes(y = ..ncount..),position = "identity",bins = 100,alpha = 0.6) + 
  facet_grid(Group ~ wavenumber,scales = "free") +
  default_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 15, angle = 45,hjust = 1,vjust = 1)) +
  xlab('Metabolic activity (CD-ratio)')

hist_CDR_plot
ggsave(filename=paste(outpath,"/","plot_hist_allwave.png", sep=""),plot=hist_CDR_plot,
       limitsize=T,width = 20,height = 20)


#####RICH

