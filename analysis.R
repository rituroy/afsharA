dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
setwd(paste(dirSrc2,"AfsharA",sep=""))


##############################################
outFormat="pdf"
outFormat="png"

##############################################
if (F) {
    clin1=read.table("docs/20171105/Book1 UM.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
    clin2=read.table("docs/Book1 UM 12-05-2017.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
    clin2=clin2[match(clin1$Patient.Number,clin2$Patient.Number),]
    k=match(names(clin1),names(clin2)); k1=which(!is.na(k)); k2=k[k1]
    for (k in 1:length(k1)) {
        if (any(is.na(clin1[,k1[k]])!=is.na(clin2[,k2[k]])) | any(clin1[,k1[k]]!=clin2[,k2[k]],na.rm=T)) {
            cat("\n\n======== ",names(clin1)[k1[k]],"\n",sep="")
            print(table(clin1=clin1[,k1[k]],clin2=clin2[,k2[k]],exclude=NULL))
        }
    }
}

if (F) {
    clin=read.table("docs/Book1 UM 12-05-2017.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
    
    clin$id=paste("X",clin$Patient.Number,sep="")
    clin$Chrom3lossBi=as.integer(clin$Chrom3loss>0)
    clin$GNAQ_GNA11=as.integer(clin$GNAQ | clin$GNA11)
    clin$BAP1_EIF1Ax_SF3B1=as.integer(clin$BAP1.loss | clin$BAP1.rearrangement | clin$EIF1Ax | clin$SF3B1)
    clin$path3=clin$Pathology; clin$path3[!clin$Pathology%in%c("1","2","3")]=NA
    clin$Chrom6pgainBi=clin$Chrom6pgain; clin$Chrom6pgainBi[!clin$Chrom6pgain%in%(0:1)]=NA
    clin$GEP3=clin$GEP; clin$GEP3[!clin$GEP%in%c("1","2","3")]=NA
    clin$tnm4=substr(clin$TNM,1,2)
}
##############################################

clin=read.table("docs/Book1 UM 122272017_deidentified.csv", sep=",", h=T, quote="", comment.char="",as.is=T,fill=T)
clin$Chrom6pgain[which(clin$CCGL.Pt.Number=="UMPA-9")]=1

#clin$id=paste("X",clin$Patient.Number,sep="")
clin$id=paste("X",1:nrow(clin),sep="")
clin$Chrom3loss=clin$Monosomy.3
clin$monosomy3Type="none"
clin$monosomy3Type[which(clin$Monosomy.3==1)]="complete"
clin$monosomy3Type[which(clin$Partial.monosomy.3==1)]="partial"
clin$monosomy3Type=clin$Monosomy.3
clin$monosomy3Type[which(clin$Partial.monosomy.3==1)]=2
clin$Chrom3lossBi=as.integer(clin$Chrom3loss>0)
clin$BAP1.loss=clin$BAP1
clin$GNAQ_GNA11=as.integer(clin$GNAQ | clin$GNA11)
clin$BAP1_EIF1Ax_SF3B1=as.integer(clin$BAP1.loss | clin$BAP1.rearrangement | clin$EIF1Ax | clin$SF3B1)
clin$path3=clin$Pathology; clin$path3[!clin$Pathology%in%c("1","2","3")]=NA
clin$Chrom6pgainBi=clin$Chrom6pgain; clin$Chrom6pgainBi[!clin$Chrom6pgain%in%(0:1)]=NA
clin$GEP3=clin$GEP; clin$GEP3[!clin$GEP%in%c("1","2","3")]=NA
clin$tnm4=substr(clin$TNM,1,2)

varInfo=data.frame(varList=c("Agem","Race","Sex","tnm4","path3","Largestbasaldiam","Thickness","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","Chrom1pLoss","Monosomy.3","Partial.monosomy.3","monosomy3Type","Chrom3loss","Chrom3lossBi","Chrom6pgainBi","Chrom6qloss","Chrom8ploss","Chrom8qgain","GNAQ","GNA11","BAP1.loss","BAP1.rearrangement","EIF1Ax","SF3B1","GNAQ_GNA11","BAP1_EIF1Ax_SF3B1","GEP3","GEP2"),
varName=c("Age","Race","Sex","TNM","Pathology","LargestBasalDiameter","Thickness","CiliaryBodyInvolvement","EOE","EpithelioidCell","Chrom1pLoss","Monosomy3","PartialMonosomy3","Chrom3Loss","Chrom3loss","Chrom3Loss","Chrom6pGain","Chrom6qLoss","Chrom8pLoss","Chrom8qGain","GNAQ","GNA11","BAP1Loss","BAP1Rearrange","EIF1Ax","SF3B1","GNAQ_GNA11","BAP1_EIF1Ax_SF3B1","GEP","GEP"),
varNameLong=c("Age","Race","Sex","TNM stage","Pathology","Largest basal tumor diameter","Tumor thickness","Ciliary body involvement","Extraocular extension","Epithelioid cells","Chrom 1p loss","Monosomy 3","Partial monosomy 3","Chrom 3 loss","Chrom 3 loss","Chrom 3 loss","Chrom 6p gain","Chrom 6q loss","Chrom 8p loss","Chrom 8q gain","GNAQ mutation","GNA11 mutation","BAP1 mutation","BAP1 rearrangement","EIF1Ax mutation","SF3B1 mutation","GNAQ_GNA11 mutation","BAP1_EIF1Ax_SF3B1 mutation","GEP","GEP"),stringsAsFactors=F)
catInfo=data.frame(catList=c("0","1"),catNameShort=c("noMut","mut"),catNameLong=c("no mutation","mutation"),stringsAsFactors=F)

fName1="_um"

save.image("tmp.RData")

##############################################

varList1=c("TNM","Pathology","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","Chrom3loss","Chrom1pLoss","Chrom6pgain","Chrom8qgain","GNAQ","GNA11","BAP1.loss","BAP1.rearrangement","EIF1Ax","SF3B1","GEP")
for (vId1 in 1:length(varList1)) {
    cat("\n\n ============ ",varList1[vId1],"\n")
    print(table(clin[,varList1[vId1]],exclude=NULL))
}

varList1=list(c("GNAQ","GNA11"),c("BAP1.loss","BAP1.rearrangement","EIF1Ax","SF3B1"))
for (vId1 in 1:length(varList1)) {
    varList2=varList1[[vId1]]
    cat("\n\n ============ ",varList2[1],"\n")
    for (vId2 in 1:(length(varList2)-1)) {
        for (vId3 in (vId2+1):length(varList2)) {
            x=table(clin[,varList2[vId2]],clin[,varList2[vId3]],exclude=NULL,dnn=c(varList2[vId2],varList2[vId3]))
            #x=table(clin[,varList2[vId2]],clin[,varList2[vId3]],exclude=NULL,dnn=varInfo$varName[match(c(varList2[vId2],varList2[vId3]),varInfo$varList)])
            cat("\n")
            print(x)
        }
    }
}

## --------------

statTbl=NULL
contingentTbl=NULL

## --------------
library(coin)

cexMain=2
cexLab=1.5
cexLab=2
varList1=c("Largestbasaldiam","Thickness")
varList2=c("Chrom3lossBi")
n=length(varList1)*length(varList2)
tmp=rep(NA,n)
tmpC=rep("",n)
tbl=data.frame(variable1=tmpC,variable2=tmpC,testType=tmpC,pv=tmp,ld=tmp,or=tmp,lcb=tmp,ucb=tmp,stringsAsFactors=F)
k=1
if (outFormat=="pdf") {
    pdf(paste("boxplot",fName1,".pdf",sep=""),width=7, height=7/2)
} else {
    png(paste("boxplot",fName1,".png",sep=""),width=2*480,height=1*480)
}
par(mfrow=c(1,2))
par(mar=c(5, 4, 4, 2) + 0.1)
par(mar=c(5, 5, 4, 2) + 0.1)
for (vId1 in 1:length(varList1)) {
    for (vId2 in 1:length(varList2)) {
        j=which(!is.na(clin[,varList1[vId1]]) & !is.na(clin[,varList2[vId2]]))
        #x=table(clin[j,varList2[vId2]])
        x=table(clin[j,varList2[vId2]],dnn=varInfo$varName[match(varList2[vId2],varInfo$varList)])
        ttl=paste(catInfo$catNameLong[match(names(x),catInfo$catList)]," (",x,")",sep="")
        #res=kruskal_test(clin[,varList1[vId1]]~as.factor(clin[,varList2[vId2]]),distribution="exact")
        #res=kruskal_test(clin[,varList1[vId1]]~as.factor(clin[,varList2[vId2]]))
        res=wilcox_test(clin[,varList1[vId1]]~as.factor(2-clin[,varList2[vId2]]),distribution="exact",conf.int=T)
        tbl$ld[k]=confint(res)$est # Hodges-Lehmann estimator
        tbl[k,c("lcb","ucb")]=confint(res)$conf.int # exact confidence interval is obtained by the algorithm described in Bauer (1972)
        tbl$variable1[k]=varList1[vId1]
        tbl$variable2[k]=varList2[vId2]
        #tbl$testType[k]="Kruskal-Wallis test"
        tbl$testType[k]="Wilcoxon test"
        tbl$pv[k]=pvalue(res)
        boxplot((clin[,varList1[vId1]]~as.factor(clin[,varList2[vId2]])),names=ttl,main=paste("Wilcoxon's test p-value ",signif(tbl$pv[k],2),sep=""),xlab=varInfo$varNameLong[match(varList2[vId2],varInfo$varList)],ylab=varInfo$varNameLong[match(varList1[vId1],varInfo$varList)],cex.main=cexMain,cex.lab=cexLab)
        k=k+1
    }
}
dev.off()
statTbl=rbind(statTbl,tbl)

median(clin[which(clin[,varList2[vId2]]==1),varList1[vId1]],na.rm=T)-median(clin[which(clin[,varList2[vId2]]==0),varList1[vId1]],na.rm=T)
res1=wilcox_test(clin[,varList1[vId1]]~as.factor(clin[,varList2[vId2]]),distribution="exact",conf.int=T)
res2=wilcox.test(clin[,varList1[vId1]]~as.factor(clin[,varList2[vId2]]),exact=T,conf.int=T)

y=1:20
x=rep(0:1,each=10); x[2:3]=1; x[14:16]=0
x=rep(0:1,each=10); x[2:3]=1; x[14:17]=0
res1=wilcox_test(y~as.factor(x),distribution="exact",conf.int=T)
res2=wilcox.test(y~as.factor(x),exact=T,conf.int=T)
res1
res2
median(y[x==1])-median(y[x==0])



## --------------

fName=paste("contingencyTbl",fName1,".txt",sep="")

varList1=c("tnm4","path3","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","Chrom1pLoss","Chrom6pgainBi","Chrom8qgain","GNAQ_GNA11","BAP1_EIF1Ax_SF3B1","GEP3")
varList1=c("tnm4","path3","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","Chrom1pLoss","Chrom6pgainBi","Chrom8qgain","GNAQ","GNA11","BAP1.loss","BAP1.rearrangement","EIF1Ax","SF3B1","GEP3")
varList2=c("Chrom3lossBi")
n=length(varList1)*length(varList2)
tmp=rep(NA,n)
tmpC=rep("",n)
tbl=data.frame(variable1=tmpC,variable2=tmpC,testType=tmpC,pv=tmp,ld=tmp,or=tmp,lcb=tmp,ucb=tmp,stringsAsFactors=F)
k=1
write.table("",file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=F)
for (vId1 in 1:length(varList1)) {
    for (vId2 in 1:length(varList2)) {
        #x=table(clin[,varList1[vId1]],clin[,varList2[vId2]])
        x=table(clin[,varList1[vId1]],clin[,varList2[vId2]],dnn=varInfo$varName[match(c(varList1[vId1],varList2[vId2]),varInfo$varList)])
        res=fisher.test(x)
        tbl$variable1[k]=varInfo$varName[match(varList1[vId1],varInfo$varList)]
        tbl$variable2[k]=varInfo$varName[match(varList2[vId2],varInfo$varList)]
        tbl$testType[k]="Fisher test"
        tbl$pv[k]=res$p.value
        if (nrow(x)==2) {
            tbl$or[k]=res$estimate
            tbl$lcb[k]=res$conf.int[1]
            tbl$ucb[k]=res$conf.int[2]
        }
        write.table(paste(c("",varInfo$varName[match(varList2[vId2],varInfo$varList)]),collapse="\t"),file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=T)
        write.table(paste(c(varInfo$varName[match(varList1[vId1],varInfo$varList)],colnames(x)),collapse="\t"),file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=T)
        write.table(x,file=fName, sep="\t", col.names=F, row.names=T, quote=F,append=T)
        write.table("",file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=T)
        k=k+1
    }
}
statTbl=rbind(statTbl,tbl)

for (k in c("pv")) statTbl[,k]=signif(statTbl[,k],4)
for (k in c("ld","or","lcb","ucb")) statTbl[,k]=round(statTbl[,k],4)
nm=c("variable1","variable2","test","p-value","location shift","odds ratio","lower confidence bound","upper confidence bound")
write.table(paste(nm,collapse="\t"),file=paste("stat",fName1,".txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=F)
write.table(statTbl,file=paste("stat",fName1,".txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=T)



## --------------
# Multivariate analysis proposed: Chromsome 3 loss vs Gene Expression Profile Class and TNM stage

varList1=c("tnm4","GEP3")
varList2=c("Chrom3lossBi")

resp=clin[,varList2]

dat=clin[,varList1]
for (vId in 1:length(varList1)) {
    model2=as.formula(paste("resp~",varList1[vId],sep=""))
    print(summary(glm(model2,family="binomial",data=dat))$coef)
}

x=clin$GEP3; x[which(!clin$GEP3%in%c("1","3"))]=NA
dat=data.frame(resp,tnm4=clin$tnm4,GEP1a2=x,stringsAsFactors=F)

x=as.integer(substr(clin$tnm4,2,2))
dat=data.frame(resp,tnmCont=x,stringsAsFactors=F)
x=clin$GEP3; x[which(!clin$GEP3%in%c("1","3"))]=NA
dat=cbind(dat,GEP1a2=x)

x=clin$tnm4; x[which(clin$tnm4%in%c("T1","T2"))]="T1T2"; x[which(clin$tnm4%in%c("T3","T4"))]="T3T4"
dat=data.frame(resp,tnm12_34=x,stringsAsFactors=F)
x=clin$GEP3; x[which(!clin$GEP3%in%c("1","3"))]=NA
dat=cbind(dat,GEP1a2=x)

x=clin$tnm4; x[which(clin$tnm4%in%c("T2"))]=NA; x[which(clin$tnm4%in%c("T3","T4"))]="T3T4"
dat=data.frame(resp,tnm1_34=x,stringsAsFactors=F)
x=clin$GEP3; x[which(!clin$GEP3%in%c("1","3"))]=NA
dat=cbind(dat,GEP1a2=x)

x=clin$tnm4; x[which(clin$tnm4%in%c("T2"))]=NA
dat=data.frame(resp,tnm134=x,stringsAsFactors=F)
x=clin$GEP3; x[which(!clin$GEP3%in%c("1","3"))]=NA
dat=cbind(dat,GEP1a2=x)

x=clin$tnm4; x[which(clin$tnm4%in%c("T2","T3"))]=NA
dat=data.frame(resp,tnm14=x,stringsAsFactors=F)
x=clin$GEP3; x[which(!clin$GEP3%in%c("1","3"))]=NA
dat=cbind(dat,GEP1a2=x)

model2=as.formula(paste("resp~tnm12_34*GEP1a2",sep=""))
model2=as.formula(paste("resp~.",sep=""))
res=glm(model2,family="binomial",data=dat)
summary(res)$coef

## ------------
## Multivariate model

library(MASS)

varListTbl=list(
c("Agem","Race","Sex","tnm4","path3","Largestbasaldiam","Thickness","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","GEP3","Chrom1pLoss","Chrom6pgainBi","Chrom6qloss","Chrom8ploss","Chrom8qgain","GNAQ_GNA11","BAP1_EIF1Ax_SF3B1"),
c("Agem","Race","Sex","tnm4","path3","Largestbasaldiam","Thickness","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","GEP3","Chrom1pLoss","Chrom6pgainBi","Chrom6qloss","Chrom8ploss","Chrom8qgain","GNAQ","GNA11","BAP1.loss","BAP1.rearrangement","EIF1Ax","SF3B1"),
c("Agem","Race","Sex","tnm4","path3","Largestbasaldiam","Thickness","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","Chrom1pLoss","Chrom6pgainBi","Chrom8qgain","GNAQ","GNA11","BAP1.loss","BAP1.rearrangement","EIF1Ax","SF3B1"),
c("tnm4","path3","Largestbasaldiam","Thickness","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","Chrom1pLoss","Chrom6pgainBi","Chrom8qgain","GNAQ_GNA11","BAP1_EIF1Ax_SF3B1","GEP3"),
c("tnm4","path3","Largestbasaldiam","Thickness","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","Chrom1pLoss","Chrom6pgainBi","Chrom8qgain","GNAQ_GNA11","BAP1_EIF1Ax_SF3B1"),
c("tnm4","path3","Largestbasaldiam","Thickness","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","Chrom1pLoss","Chrom6pgainBi","Chrom8qgain","GNAQ","GNA11","BAP1.loss","BAP1.rearrangement","EIF1Ax","SF3B1"),
c("tnm4","path3","Largestbasaldiam","Thickness","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","Chrom1pLoss","Chrom6pgainBi","Chrom8qgain","GNAQ","GNA11","BAP1.loss","BAP1.rearrangement","EIF1Ax","SF3B1","GEP3")
)
varList2=c("Chrom3lossBi")

for (k in 1:length(varListTbl)) {
    varListTbl[[k]]=varListTbl[[k]][which(!varListTbl[[k]]%in%c("BAP1.rearrangement"))]
}

for (modelFlag in c("gep1a2","gep1a1b2")) {
    fName=paste("multivariate_",modelFlag,fName1,".txt",sep="")
    write.table("Multivariate model with Chrom3Loss as the response variable",file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=F)
    for (vId in 1:length(varListTbl)) {
        varList1=varListTbl[[vId]]
        dat=clin[,c(varList2,varList1)]
        if (modelFlag=="gep1a2") {
            if (length(grep("GEP",varList1))==0) next
            varList1=sub("GEP3","GEP2",varList1)
            dat$GEP2=clin$GEP3; dat$GEP2[which(dat$GEP2=="2")]=NA
            dat=dat[,c(varList2,varList1)]
        }

        write.table(paste("\nVariables considered: ",paste(varInfo$varName[match(varList1,varInfo$varList)],collapse=", "),sep=""),file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=T)
        j=1:nrow(dat)
        for (k in 1:ncol(dat)) {
            j=j[!is.na(dat[j,k])]
        }
        dat=dat[j,]
        #model2=as.formula(paste(varList2,"~.",sep=""))
        names(dat)=varInfo$varName[match(names(dat),varInfo$varList)]
        model2=as.formula(paste(varInfo$varName[match(varList2,varInfo$varList)],"~.",sep=""))

        fit=glm(model2,family="binomial",data=dat)
        fit2=stepAIC(fit,k=log(nrow(dat)))
        res=summary(fit2)
        write.table(paste("Final model: ",varList2," ~ ",fit2$formula[3],sep=""),file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=T)
        tbl=as.data.frame(res$coef)
        for (k in c("Estimate","Std. Error","z value")) tbl[,k]=signif(tbl[,k],2)
        for (k in "Pr(>|z|)") tbl[,k]=signif(tbl[,k],2)
        tbl=cbind(rownames(res$coef),tbl)
        names(tbl)=c("",colnames(res$coef))
        write.table(tbl,file=fName, sep="\t", col.names=T, row.names=F, quote=F,append=T)
    }
}

##############################################
##############################################
## NOT USED

if (F) {
    library(foreign)

    tbl1=read.dta("PB_102517.dta",convert.factors=F)

    library(readstata13)

    dat <- read.dta13(system.file("extdata/statacar.dta", package="readstata13"))

    tbl2=read.dta("armin_103017.dta",convert.factors=F)
    tbl3=read.dta("armin_time_103017.dta",convert.factors=F)
}





##############################################
