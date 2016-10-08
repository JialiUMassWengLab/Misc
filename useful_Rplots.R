plotPhysicianRating <- function() {
data <- read.csv("Non-Opioid-physician-rating.csv",header=FALSE,strip.white = TRUE, stringsAsFactors = FALSE)
names(data) <- c("pid","clinic","visit","Explain inadequacy","Explain side-effects",
                 "Change medication","Change dosage","Reduce trial-and-error","rate")

pdf("physician_rating_pie.pdf")
par(mar=c(2,2,5,2), cex=1.3)
pie(table(data[data$visit==2,]$rate),edges=3600,radius=0.7)
legend("topright",c("1: Not at all valuable","2: Slightly valuable",
                    "3: Moderately valuable","4: Mostly valuable","5: Extremely valuable"),bty="n",cex=0.65)
dev.off()
}

plotPhysicianAnswer <- function() {
data <- read.csv("Non-Opioid-physician-rating.csv",header=FALSE,strip.white = TRUE, stringsAsFactors = FALSE)
names(data) <- c("pid","clinic","visit","Explain inadequacy","Explain side-effects",
                 "Change medication","Change dosage","Reduce trial-and-error","rate")

data2 <- data[data$visit==2,]
matrix <- rbind(table(data2[data2[,4]!="" & data2[,4]!="N/A (test not ordered)",4]),table(data2[data2[,5]!="" & data2[,5]!="N/A (test not ordered)",5]),table(data2[data2[,6]!="" & data2[,6]!="N/A (test not ordered)",6]),
                table(data2[data2[,7]!="" & data2[,7]!="N/A (test not ordered)",7]),table(data2[data2[,8]!="" & data2[,8]!="N/A (test not ordered)",8]))
print(matrix)
matrix <- matrix/apply(matrix,1,sum)
pdf("physician_answers_bar.pdf", 21, 0.8*nrow(matrix)+6)
par(mar=c(3,15,3,2)+0.1, fin=c(21,0.8*nrow(matrix)+6), pin=c(10,0.8*nrow(matrix)), cex.axis=1.5, cex.lab=1.5)
a <- barplot(t(matrix),horiz=TRUE,width=0.5,space=0.8,yaxt="n", ylim=c(0,4), xlab="Fraction")
text(cex=1.5, font=2, x=-0.01*(par("usr")[2]-par("usr")[1]), y=a, names(data2)[4:8], xpd=TRUE, pos=2)
#legend(0.9,4.5,c("No","Yes"),fill=c("black","grey"))
dev.off()
}


plotMED_vs_SNP <- function(snp,drug_list) {
  source("Non-Opioid-Response.R")
  snps <- getSNP("all_SNP_calls.csv")
  resp <- getMED2(drug_list)
  types <- read.csv("reference_files/Allele_class.csv")
  mat <- joinGenotypeMED(snps[snps$rs==snp,],resp,types[[snp]][2])
  major_allele <- levels(factor(mat[mat$allele_type==FALSE,]$call))[1]
  minor_allele <- levels(factor(mat[mat$allele_type==TRUE,]$call))
  minor_allele <- paste(minor_allele[1],minor_allele[2],sep=" or ")
  
  matrix <- table(mat$allele_type, mat$grade)
  print(matrix)
  matrix <- matrix/apply(matrix,1,sum)
  t <- matrix[,1]
  matrix[,1] <- matrix[,2]
  matrix[,2] <- t
  
  main_t <- paste("MED response to Gabapentin",snp,sep="\n")
  pdf(paste(snp,"_MED_response.pdf",sep=""), 18, 7)
  par(mar=c(3,12,3,2)+0.1, fin=c(18,7), pin=c(10,2), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
  barplot(t(matrix),horiz=TRUE,las=1,names.arg=c(major_allele,minor_allele),width=0.5,space=1,ylim=c(0,2),xlab="Fraction",main=main_t)
  dev.off()
}

plotReportMED <- function(drug) {
  source("Non-Opioid-Response.R")
  m <- joinReportMED(drug)
  matrix <- table(m$isGoodResponder,m$grade)
  matrix <- matrix/apply(matrix,1,sum)
  t <- matrix[,1]
  matrix[,1] <- matrix[,2]
  matrix[,2] <- t
  
  main_t <- "MED response to Gabapentin VS Proove PBIO13 test"
  pdf("Gabapentin_MED_report.pdf", 18, 7)
  par(mar=c(3,12,3,2)+0.1, fin=c(18,7), pin=c(10,2), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
  barplot(t(matrix),horiz=TRUE,las=1,names.arg=c("Bad responder","Good responder"),width=0.5,space=1,ylim=c(0,2),xlab="Fraction",main=main_t)
  dev.off()
}

plotPhysicianImprovement <- function() {
  library(VennDiagram)
  pdf("Intervention_Venn.pdf")
  grid.newpage()
  grid.text("n=409",0.95,0.9,gp=gpar(fontsize=15,fontface="bold"))
  #grid.rect(width=0.98, height=0.98, gp=gpar(fill=NA,lwd=2))
  draw.pairwise.venn(31,97,12,scaled=FALSE,category=c("Physician Intervened","Patient Response Improved"),
                     lty=rep("blank",2),fill=c("light blue","pink"), alpha=rep(0.5,2),cat.pos=c(15,-15),
                     fontface=rep("bold",3),fontsize=rep(20,3),cat.fontface=rep("bold",2),cat.fontsize=rep(20,2))
  dev.off()
}

plotLegend <- function() {
  pdf("legend.pdf")
  plot(c(1:10),c(1:10),type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  legend("topright",c("No","Yes"),fill=c("black","light grey"),text.font=2,cex=1.5)
  legend("topleft",c("Unguided","Guided"),fill=c("black","light grey"),text.font=2,cex=1.5)
  legend("bottomleft",c("Good","Average","Poor"),fill=c("black","dark grey","light grey"),text.font=2,cex=1.5)
  dev.off()
}

plotGuidedVsPhysicianRating <- function() {
    source("Non-Opioid-Response.R")
    ratings <- getPhysicianRating()
    ratings <- ratings[ratings$visit==3,]
    print(wilcox.test(ratings[ratings$guided=="Yes","rate"],ratings[ratings$guided=="No","rate"]))
    matrix <- table(ratings$guided,ratings$rate)
    print(matrix)
    matrix <- matrix/apply(matrix,1,sum)
    pdf("Guided_vs_PhysicianRating_bar.pdf")
    barplot(matrix,beside=TRUE,xlab="Physician Rating",ylab="Fraction",ylim=c(0,0.5))
    dev.off()
}

#plotPhysicianAnswer()
#plotMED_vs_SNP("rs242939",gabapentin)
#plotMED_vs_SNP("rs7016778",gabapentin)
#plotReportMED(gabapentin)
#plotPhysicianImprovement()
#plotLegend()
plotGuidedVsPhysicianRating()
