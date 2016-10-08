getMED <- function () {
    data<-read.csv("all_r_m_aug1.MED.csv",header=FALSE,strip.white = TRUE, stringsAsFactors = FALSE)
    names(data)<-c("pid","visit","drug","p.rate")

    resp2<-data[data$visit==2,]
    resp2["pair"]<-paste(resp2$pid,resp2$drug,sep=".")
    resp3<-data[data$visit==3,]
    resp3["pair"]<-paste(resp3$pid,resp3$drug,sep=".")
    
    resp<-merge(resp2,resp3,by="pair")
    resp<-resp[,c("pair","pid.x","drug.x","p.rate.x","p.rate.y")]
    names(resp)<-c("pair","pid","drug","p.rate.visit2","p.rate.visit3")
    
    return(resp)
}

getMED2 <- function(drug_list) {
    MED<-read.csv("all_r_m_aug1.MED.csv",header=FALSE,strip.white = TRUE, stringsAsFactors = FALSE)
    names(MED)<-c("pid","visit","drug","p.rate")
    MED <- MED[!is.na(MED$p.rate) & MED$p.rate!=""
               & MED$p.rate!="na",]

    MED <- data.frame(do.call(rbind,apply(MED,1,function (x) {
        for (y in drug_list) {
            if (grepl(y,x[3])) {return(as.list(x))}
        }
    })))    
    
    MED$p.rate <- as.numeric(MED$p.rate)
    MED$pid <- as.character(MED$pid)
    MED <- MED[!is.na(MED$p.rate) & MED$p.rate>=0 & MED$p.rate<=5,]
    MED <- MED[with(MED, order(pid,-p.rate)),]

    resp2 <- unique(MED[MED$visit==2,c("pid","p.rate")])
    resp2 <- resp2[!duplicated(resp2$pid),]
    #resp2 <- resp2[!resp2$pid %in% resp2$pid[which(duplicated(resp2$pid))],]
    resp1 <- unique(MED[!(MED$pid %in% resp2$pid) & MED$visit==1,c("pid","p.rate")])
    resp1 <- resp1[!duplicated(resp1$pid),]
    #resp1 <- resp1[!resp1$pid %in% resp1$pid[which(duplicated(resp1$pid))],]
    resp3 <- unique(MED[!(MED$pid %in% resp2$pid) & !(MED$pid %in% resp1$pid)
                 & MED$visit==3,c("pid","p.rate")])
    resp3 <- resp3[!duplicated(resp3$pid),]
    #resp3 <- resp3[!resp3$pid %in% resp3$pid[which(duplicated(resp3$pid))],]
    resp4 <- unique(MED[!(MED$pid %in% resp2$pid) & !(MED$pid %in% resp1$pid)
                 & !(MED$pid %in% resp3$pid) & MED$visit==4,c("pid","p.rate")])
    resp4 <- resp4[!duplicated(resp4$pid),]
    #resp4 <- resp4[!resp4$pid %in% resp4$pid[which(duplicated(resp4$pid))],]
    
    resp <- rbind(resp1,resp2,resp3,resp4, stringAsFactors=FALSE, make.row.names=FALSE)

    return(resp)
}

getNRS <- function() {
    data <- read.csv("all_rm_sep6.NRS.csv",header=FALSE,strip.white = TRUE, stringsAsFactors = FALSE)
    names(data) <- c("pid","visit","before","after")
    #print(dim(data))
    data <- data[data$before > 0,]
    data["change"] <- as.vector(apply(data, 1, function(x) {
        return(as.numeric(x[3])-as.numeric(x[4]))}))
    data <- data[,c("pid","visit","before")]
    #print(dim(data))
    return(data)
}

getPhysicianRating <- function() {
    data <- read.csv("Non-Opioid-physician-rating.csv",header=FALSE,strip.white = TRUE, stringsAsFactors = FALSE)
    names(data) <- c("pid","clinic","visit","Explain inadequacy","Explain side-effects",
                                      "Change medication","Change dosage","Reduce trial-and-error","rate")
    data["guided"] <- apply(data,1,function(x) {
        if (x[6]=="Yes" | x[7] == "Yes" | x[8] == "Yes") {return("Yes")}
        else {return("No")}
    })
    return(data)
}

getSNP <- function (fn) {
    snps <- read.csv(fn,header=FALSE,stringsAsFactors = FALSE)
    names(snps) <- c("rs","call","sid","pid")
    snps <- snps[!grepl("\\d",snps$call)
                 & snps$call!="" & snps$call!="NOAMP"
                 & snps$call!="INV" & snps$call!="UND",]
    return(snps)
}

getPatientInfo <- function() {
    demo <- read.csv("reference_files/patientInfo.csv",header=FALSE,stringsAsFactors=FALSE)
    names(demo) <- c("pid","gender","race")
    demo <- unique(demo)
    demo <- demo[!demo$pid %in% demo$pid[which(duplicated(demo$pid))],]
    return(demo)
}

joinGenotypeMED <- function(genotype, MED, minor) {
    genotype <- unique(genotype[,c("pid","call")])
    genotype <- genotype[!genotype$pid %in% genotype$pid[which(duplicated(genotype$pid))],]
    genotype["allele_type"] <- as.vector(apply(genotype,1,function(x) {grepl(minor,x[2])}))

    m <- merge(genotype,MED,by="pid")
    m["grade"] <- as.vector(apply(m,1,function(x) {
        if (x[4] < 2) {"Poor"}
        else if (x[4] < 4) {"Average"}
        else {"Good"}
    }))
    return(m)
}

testResponseVsPhysician <- function(drug_list) {
    data<-read.csv("Non-Opioid-eval.csv",header=FALSE,strip.white = TRUE, stringsAsFactors = FALSE)
    names(data)<-c("pid","visit","drug","change","d.rate")
    ch2<-data[data$visit==2,c("pid","drug","change")]
    ch2["pair"]<-paste(ch2$pid,ch2$drug,sep=".")
    rate3<-data[data$visit==3,c("pid","drug","d.rate")]
    rate3["pair"]<-paste(rate3$pid,rate3$drug,sep=".")
    physician<-merge(ch2,rate3,by="pair")
    physician<-physician[,c("pair","pid.x","drug.x","change","d.rate")]
    names(physician)<-c("pair","pid","drug","change","d.rate")
    physician<-physician[physician$change!="" & physician$d.rate!="" 
                         & !is.na(physician$change) & !is.na(physician$d.rate),]

    resp <- getMED()
    #result<-merge(resp,physician,by="pair")
    result<-merge(resp,ch2,by="pair")
    result<-result[,!names(result) %in% c("pid.y","drug.y")]
    if (length(drug_list)>1 | drug_list != "ALL") {
        result <- result[result$drug.x %in% drug_list,]
    }
    result$p.rate.visit2<-as.numeric(result$p.rate.visit2)
    result$p.rate.visit3<-as.numeric(result$p.rate.visit3)
    result<-result[!is.na(result$p.rate.visit2) & !is.na(result$p.rate.visit3),]
    
    result["intervention"]<-!grepl("No Change",result$change)
    result["p.improve"]<-(result$p.rate.visit3 > result$p.rate.visit2)
    ctable <- table(result$intervention,result$p.improve)
    print(ctable)
    phyper(ctable[4],ctable[3]+ctable[4],ctable[1]+ctable[2],ctable[2]+ctable[4],lower.tail=FALSE)
}

printDetails <- function(snp,drug_list) {
    snps <- getSNP("all_SNP_calls.csv")
    resp <- getMED2(drug_list)
    types <- read.csv("reference_files/Allele_class.csv")
                
    mat <- joinGenotypeMED(snps[snps$rs==snp,],resp,types[[snp]][2])
    print(dim(mat))
    print(chisq.test(mat$allele_type, mat$grade))
    print(table(mat$allele_type, mat$grade))
    print(table(mat$call, mat$grade))
    print(table(mat$allele_type, as.numeric(mat$p.rate)))
    print(table(mat$call, as.numeric(mat$p.rate)))
    
    demo <- getPatientInfo()
    mat1 <- merge(demo,mat,by="pid")
    print(dim(mat1))
    #print(table(mat1$gender,mat1$grade))
    #print(chisq.test(mat1$gender,mat1$grade))
    #mat2 <- mat1[mat1$race %in% c("W","B/AA","A","H/L"),]
    #print(table(mat2$race,mat2$grade))
    #print(chisq.test(mat2$race,mat2$grade))
    #mat3 <- mat1[mat1$race=="H/L",]
    #print(table(mat3$allele_type,mat3$grade))
    #print(chisq.test(mat3$allele_type,mat3$grade))

    test <- multinom(grade ~ allele_type + gender + race, data=mat1)
    print(summary(test))
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p <- (1 - pnorm(abs(z), 0, 1))*2
    print(p)
    print(min(p[3:4]))
}

testGeneticAssoc <- function(drug_list) {
    snps <- getSNP("all_SNP_calls.csv")
    MED <- getMED2(drug_list)
    demo <- getPatientInfo()
    
    types <- read.csv("reference_files/Allele_class.csv")    
    all_snps <- unique(names(types))
    #info <- read.csv("reference_files/SNPinfo.csv",header=FALSE)
    #nonCYP <- info[!grepl("CYP",info$V2),1]
    filter <- read.table("reference_files/SNPs_used_PBIO2_PBIO4",header=FALSE)
    nonCYP <- filter$V1
    
    pv <- as.numeric()
    for (x in all_snps) {
        if (!x %in% nonCYP) {
            print(x)
            merged <- joinGenotypeMED(snps[snps$rs==x,],MED,types[[x]][2])
            merged <- merge(demo,merged,by="pid")
            print(dim(merged))
            if (length(levels(factor(merged$allele_type))) >= 2
                #& length(levels(factor(merged[merged$allele_type==TRUE,]$grade))) == 3
                #& length(levels(factor(merged[merged$allele_type==FALSE,]$grade))) == 3
                ) {
                #chisq_pv <- chisq.test(merged$allele_type,merged$grade)$p.value
                #print(chisq_pv)
                #pv <- append(pv,chisq_pv)
                test <- multinom(grade ~ allele_type + gender + race, data=merged)
                z <- summary(test)$coefficients/summary(test)$standard.errors
                p <- (1 - pnorm(abs(z), 0, 1))*2
                logistic_pv <- min(p[3:4])
                print(logistic_pv)
                pv <- append(pv,logistic_pv)
            }
        }
    }
    print(pv)
    print(p.adjust(pv,method="fdr"))
}

joinReportMED <- function(drug) {
    MED <- getMED2(drug)
    demo <- getPatientInfo()
    filename <- paste(tolower(drug[1]),"Response.csv",sep="")
    print(filename)
    report <- read.csv(filename, header=FALSE, stringsAsFactors=FALSE)
    names(report) <- c("sid","pid","rank","maxRank","isGoodResponder")
    report$rank <- as.numeric(report$rank)
    report <- report[!is.na(report$rank),]
    report <- unique(report[with(report, order(pid,-rank)),])
    report <- report[!duplicated(report$pid),c("pid","rank","isGoodResponder")]
    m <- merge(report, MED, by="pid")
    m["grade"] <- as.vector(apply(m,1,function(x) {
        if (x[4] < 2) {"Poor"}
        else if (x[4] < 4) {"Average"}
        else {"Good"}
    }))
    m <- merge(m, demo, by="pid")
    return(m)
}

testReportMED <- function(drug) {
    m <- joinReportMED(drug)
    print(dim(m))

    print(table(m$isGoodResponder,m$grade))
    print(chisq.test(m$isGoodResponder,m$grade))
    #print(cor.test(m$p.rate,m$rank,method="spearman"))
    #print(wilcox.test(m[m$isGoodResponder=="True","p.rate"],m[!m$isGoodResponder=="False","p.rate"]))

    test <- multinom(grade ~ isGoodResponder + gender + race, data=m)
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p <- (1 - pnorm(abs(z), 0, 1))*2
    print(summary(test))
    print(p)
    print(min(p[3:4]))
}

testReportMEDall <- function() {
    m <- joinReportMED(gabapentin)
    m <- rbind(m, joinReportMED(ibuprofen))
    m <- rbind(m, joinReportMED(duloxetine))
    m <- rbind(m, joinReportMED(acetaminophen))
    m <- rbind(m, joinReportMED(alprazolam))
    print(dim(m))

    print(table(m$isGoodResponder,m$grade))
    print(chisq.test(m$isGoodResponder,m$grade))
    #print(cor.test(m$p.rate,m$rank,method="spearman"))
    #print(wilcox.test(m[m$isGoodResponder=="True","p.rate"],m[!m$isGoodResponder=="False","p.rate"]))

    test <- multinom(grade ~ isGoodResponder + gender + race, data=m)
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p <- (1 - pnorm(abs(z), 0, 1))*2
    print(summary(test))
    print(p)
    print(min(p[3:4]))
}

testGuidedVsPatientNRS <- function() {
    ratings <- getPhysicianRating()
    pain <- getNRS()
    guided <- vector()
    for (x in c(1:nrow(ratings))) {if (ratings[x,"guided"]=="Yes" & !ratings[x,"pid"] %in% guided) guided <- append(guided,ratings[x,"pid"])}
    after <- pain[pain$visit==3,]
    before <- pain[pain$visit==2,]
    before <- rbind(before, pain[!pain$pid %in% before$pid & pain$visit==1,])
    merged <- merge(before,after,by="pid")
    merged["effect"] <- merged[,5]-merged[,3]
    guided_p <- merged[merged$"pid" %in% guided,]
    nonguided_p <- merged[!merged$"pid" %in% guided,]

    print(wilcox.test(guided_p$effect,nonguided_p$effect))    
}



nonOpioidDrugs <- c("IBUPROFEN","GABAPENTIN","DULOXETINE","ACETAMINOPHEN","ALPRAZOLAM","CYMBALTA",
                    "TYLENOL","PANADOL","ADVIL","BRUFEN","MORTIN","NEROFEN","FANATREX","GABARONE",
                    "GRALISE","NEUROTIN","NUPENTIN","NEOGAB","XANAX")
opioidDrugs <- c("OXYCODONE","HYDROCODONE","TRAMADOL","MORPHINE","HYDROMORPHONE","ROXICODONE",
                 "OXYCOTIN","OXECTA","OXYIR","ENDONE","OXYNORM","XODOL","VICODIN","NORCO",
                 "LORTAB","LORCET","RYZOLT","ULTRAM","CONZIP","AVINZA","DURAMORPH","KADIAN",
                 "DEPODUR","ASTRAMORPH","EXALGO","DILAUDID","DILAUDID HP")
gabapentin <- c("GABAPENTIN","FANATREX","GABARONE","GRALISE","NEUROTIN","NUPENTIN","NEOGAB")
ibuprofen <- c("IBUPROFEN","ADVIL","BRUFEN","MOTRIN","NEUROFEN")
acetaminophen <- c("ACETAMINOPHEN","TYLENOL","PANADOL")
alprazolam <- c("ALPRAZOLAM","XANAX")
duloxetine <- c("DULOXETINE","CYMBALTA")

library(nnet)

#testResponseVsPhysician(nonOpioidDrugs)
#testResponseVsPhysician(opioidDrugs)
#testResponseVsPhysician("ALL")

#testGeneticAssoc(c("IBUPROFEN","ADVIL","BRUFEN","MOTRIN","NEUROFEN"))
#testGeneticAssoc(c("GABAPENTIN","FANATREX","GABARONE","GRALISE","NEUROTIN","NUPENTIN","NEOGAB"))
#testGeneticAssoc(c("ACETAMINOPHEN","TYLENOL","PANADOL"))
#testGeneticAssoc(c("ALPRAZOLAM","XANAX"))
#testGeneticAssoc(c("DULOXETINE","CYMBALTA"))

#printDetails("rs242939",gabapentin)
#printDetails("rs7016778",gabapentin)
#printDetails("rs28365083",c("DULOXETINE","CYMBALTA"))

#testReportMED(gabapentin)
#testReportMED(ibuprofen)
#testReportMED(acetaminophen)
#testReportMED(alprazolam)
#testReportMED(duloxetine)
#testReportMEDall()
