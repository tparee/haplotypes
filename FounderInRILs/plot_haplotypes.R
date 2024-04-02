##################################
### Plot haplotype

library(ggplot2)
genomatrix = read_tsv("~/Desktop/happlotype_TUC/Genotype_Matrix.tsv")
genomatrix = as.data.frame(genomatrix)
genomatrix = subset(genomatrix, CHROM == "I")
snps = genomatrix[,1:3]
colnames(snps) = tolower(colnames(snps))
genomatrix = genomatrix[,5:ncol(genomatrix)]
colnames(genomatrix)

foundernames = c("AB1","CB4507","CB4852","CB4855","CB4856","CB4858","JU319","JU345","JU400","MY1","MY16","N2","PB306","PX174","PX179","RC301")

# Change in -1/0/1 in 0/0.5/1 (I think it should wroks bot way but not sure)
allnames = colnames(genomatrix)
genomatrix = do.call(cbind,lapply(1:ncol(genomatrix), function(i){as.numeric(genomatrix[,i])}))
genomatrix[genomatrix==0] = 0.5
genomatrix[genomatrix==-1] = 0
colnames(genomatrix)=allnames

founders = genomatrix[,colnames(genomatrix) %in% foundernames]
rils = genomatrix[,!(colnames(genomatrix) %in% foundernames)]


########################################################################
############# Plot mosaic of haplotypes ################################
########################################################################

haplop = do.call(rbind,lapply(1:nrow(haplotype), function(h){
  founder = paste(LETTERS[which(haplotype[h,3:(ncol(haplotype)-1)] > 0)], collapse = ";")
  data.frame(pos1=haplotype[h,1], pos2 = haplotype[h,2], founder = founder, rilname = haplotype[h,ncol(haplotype)])
}))


#haplocolors = c("black","#C0274C", "#C92A1D", "#F09F47", "#E74A26", "#E3E14A", "#F8DE7A", "#90C98C", 
#         "#BEE5ED", "#94C0E4", "#6A9ECA", "#3D57A6", "#3754A0", "#302275", "#32308C", "#3C8885","#3A9378")[1:ncol(founders)]

haplocolors = c("black","#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF",
                       "#E18726", "#332F92", "#FFDC91", "#008886", "#90C98C")
                       


colortoplot = as.character(factor(haplop$founder, levels = LETTERS[1:ncol(founders)], labels = haplocolors))
#colortoplot[is.na(colortoplot)] = "#B09C85FF"#"#f6f8e3"

haplop$haplocolor  = colortoplot


haplop$y = as.numeric(factor(haplop$rilname))

library(ggplot2)

ggplot(haplop)+theme_base()+
  geom_rect(aes(xmin=pos1/1e6,xmax=pos2/1e6, ymin=y, ymax=y+0.90,fill=haplocolor), color = NA)+
  scale_fill_manual(values = unique(haplop$haplocolor), breaks = unique(haplop$haplocolor))+
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0,0), breaks = seq(1,16,1), lim = c(-0.2, 15.2))+
  scale_y_continuous(expand = c(0,0), breaks = c(27,77)+7, labels = c("Mutant","Wild-type"),lim = c(0, max(haplop$y)+9))+
  theme(axis.text.y = element_text(angle = 90, size = 17),
        axis.line.y = element_line(color = "white"),
        axis.ticks.y = element_line(color = "white"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15))+ylab("")+xlab("")+
  geom_segment(x=-0.1,xend=-0.1,y=0.5,yend=50.5, linewidth = 0.7)+
  geom_segment(x=-0.1,xend=-0.1,y=52.5,yend=102, linewidth = 0.7)+
  xlab("Physical size (Mb)")+
  coord_cartesian(xlim = c(1,3))


