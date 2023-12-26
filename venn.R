#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af

library(randomcoloR)
library(venn) 

#文件名
A="hubGenes"
B="IL-17 signaling pathway"
C="Rheumatoid arthritis"
D="NF-kappa B"
E="Bladder cancer"
G="KSHV"

#构建一个列表
geneList=list()
rt=read.table(paste0(A,".txt"),header=F,sep="\t",check.names=F)      
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)                 
geneList[[A]]=uniqGene                   
uniqLength=length(uniqGene)
print(paste("1",uniqLength,sep=" "))
rt=read.table(paste0(B,".txt"),header=F,sep="\t",check.names=F)    
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[B]]=uniqGene
uniqLength=length(uniqGene)
print(paste("2",uniqLength,sep=" "))
rt=read.table(paste0(C,".txt"),header=F,sep="\t",check.names=F)    
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[C]]=uniqGene
uniqLength=length(uniqGene)
print(paste("3",uniqLength,sep=" "))
rt=read.table(paste0(D,".txt"),header=F,sep="\t",check.names=F)      
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)                 
geneList[[D]]=uniqGene                   
uniqLength=length(uniqGene)
print(paste("4",uniqLength,sep=" "))

rt=read.table(paste0(E,".txt"),header=F,sep="\t",check.names=F)      
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)                 
geneList[[E]]=uniqGene                   
uniqLength=length(uniqGene)
print(paste("5",uniqLength,sep=" "))

rt=read.table(paste0(G,".txt"),header=F,sep="\t",check.names=F)      
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)                 
geneList[[G]]=uniqGene                   
uniqLength=length(uniqGene)
print(paste("6",uniqLength,sep=" "))

mycol <- distinctColorPalette(6)
pdf(file="hub.pdf",width=5,height=5)                                                
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
dev.off()

intersectGenes=Reduce(intersect,geneList)          
write.table(file="hub.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F) 

#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
