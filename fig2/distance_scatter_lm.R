Args <- commandArgs() ##R --no-save < cmd.r otu_table.xls env.txt
#rm (list=(ls(all=T)))
#GuoqiLiu
#liuguoqii@163.com
# https://cloud.tencent.com/developer/article/1666820
# 测试数据的来源
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4845060/
library("vegan")##This is vegan 2.5-5
library("ggplot2")#‘1.5.10’
library("ggplot2")
library(ggplot2)
library(ggpubr)
library(basicTrendline)
####
#计算环境因子的欧氏距离
site <- read.delim(Args[4], sep = '\t', row.names = 1, check.names = FALSE)
#dist(x, method = "euclidean"
site_dis  <- as.matrix(dist(scale(site), method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
site_dis_D<-as.data.frame(site_dis)
site_dis_D$`ID` <- rownames(site_dis_D)
site_dis_D <- site_dis_D[,c("ID",colnames(site_dis_D)[1:(ncol(site_dis_D)-1)])]
write.table(site_dis_D ,file="euclidean.txt",sep="\t",quote=F,row.names=F)

site_dis[upper.tri(site_dis)] <- 0
#将采样点地理距离矩阵转换为两两对应数值的数据框结构
site_dis <- reshape2::melt(site_dis)
site_dis <- subset(site_dis, value != 0)
head(site_dis)

###########vegdist里存放的是表示群落之间不相似性程度的值，典型的是Bray-Curtis dissimilarity。#######
spe <- read.delim(Args[3], sep = '\t', row.names = 1, check.names = FALSE)
spe <- data.frame(t(spe))
#vegan 包 vegdist() 计算群落间物种组成 Bray-curtis 相异度矩阵

comm_sim <- as.matrix(vegan::vegdist(spe, method = 'bray'))
comm_sim_D <- as.data.frame(comm_sim )
comm_sim_D$`ID` <- rownames(comm_sim_D)
comm_sim_D <- comm_sim_D[,c("ID",colnames(comm_sim_D)[1:(ncol(comm_sim_D)-1)])]

write.table(comm_sim_D ,file="bray.txt",sep="\t",quote=F,row.names=F)

comm_sim2 <- as.matrix(vegan::vegdist(spe, method = 'jaccard'))
comm_sim2_D <- as.data.frame(comm_sim2 )
comm_sim2_D$`ID` <- rownames(comm_sim2_D)
comm_sim2_D <- comm_sim2_D[,c("ID",colnames(comm_sim2_D)[1:(ncol(comm_sim2_D)-1)])]
write.table(comm_sim2_D ,file="jaccard.txt",sep="\t",quote=F,row.names=F)
#将矩阵转换为两两群落对应数值的数据框结构
diag(comm_sim) <- 0  #去除群落相似度矩阵中的对角线值，它们是样本的自相似度
comm_sim[upper.tri(comm_sim)] <- 0  #群落相似度矩阵是对称的，因此只选择半三角（如下三角）区域的数值即可
comm_sim <- reshape2::melt(comm_sim)
comm_sim <- subset(comm_sim, value != 0)
head(comm_sim)
###########################################
diag(comm_sim2) <- 0  #去除群落相似度矩阵中的对角线值，它们是样本的自相似度
comm_sim2[upper.tri(comm_sim2)] <- 0  #群落相似度矩阵是对称的，因此只选择半三角（如下三角）区域的数值即可
comm_sim2 <- reshape2::melt(comm_sim2)
comm_sim2 <- subset(comm_sim2, value != 0)
head(comm_sim2)
#采样点距离和群落相似度数据合并
comm_dis <- merge(comm_sim, site_dis, by = c('Var1', 'Var2'))
names(comm_dis) <- c('site1', 'site2', 'comm_sim', 'site_dis')
###########################################
comm_dis2 <- merge(comm_sim2, site_dis, by = c('Var1', 'Var2'))
names(comm_dis2) <- c('site1', 'site2', 'comm_sim', 'site_dis')
write.table(comm_dis,file="comm_dis_Bray-Curtis.txt",sep="\t",quote=F,row.names=F)
write.table(comm_dis2,file="comm_dis_Jaccard.txt",sep="\t",quote=F,row.names=F)
mycol <- c(34, 51, 142, 26, 31, 371, 36, 7, 12, 30, 84, 88, 116, 121, 77, 56, 386, 373, 423, 435, 438, 471, 512, 130, 52, 47, 6, 11, 43, 54, 367,
 382, 422, 4, 8, 375, 124, 448, 419, 614, 401, 403, 613, 583, 652, 628, 633, 496, 638, 655, 132, 503, 24)
mycol <-colors()[rep(mycol,30)]
pdf (paste(Args[3],"_Bray_Curtis.scatter.pdf",sep=""),width=6,height=6)
#plot(0:max(comm_dis$site_dis),0:max(comm_dis$comm_sim),xaxt="n",yaxt="n",type="n",xlab="",ylab="",xaxs="i",yaxs="i")
#abline(0,1)
 #axis(2,at=seq(0,1,by=.2))
 #axis(1,at=seq(0,1,by=.2))
par(new=T)
trendline(comm_dis$site_dis, comm_dis$comm_sim,model="line2P", ePos.x=4,ePos.y=.2,
xaxt="n",yaxt="n",summary=TRUE, eDigit=5,show.equation=FALSE,
linecolor="red",lwd=4,lty=1,pch=16,cex=0.6,xlab=("euclidean dissimilarity"),ylab=("Bray-Curtis dissimilarity"))
axis(2);axis(1);
#abline(lm(lm(y2018~y2017,data=da.all)),lwd=1,lty=2,cex=0.6)
dev.off()
#########################################
pdf (paste(Args[3],"_Jaccard.scatter.pdf",sep=""),width=6,height=6)
#plot(0:max(comm_dis2$site_dis),0:max(comm_dis2$comm_sim),xaxt="n",yaxt="n",type="n",xlab="",ylab="",xaxs="i",yaxs="i")
#abline(0,1)
 #axis(2,at=seq(0,1,by=.2))
 #axis(1,at=seq(0,1,by=.2))
par(new=T)
trendline(comm_dis2$site_dis, comm_dis2$comm_sim,model="line2P", ePos.x=4,ePos.y=.2,
xaxt="n",yaxt="n",summary=TRUE, eDigit=5,show.equation=FALSE,
linecolor="red",lwd=4,lty=1,pch=16,cex=0.6,xlab=("euclidean dissimilarity"),ylab=("Jaccard dissimilarity"))
axis(2);axis(1);
#abline(lm(lm(y2018~y2017,data=da.all)),lwd=1,lty=2,cex=0.6)
dev.off()