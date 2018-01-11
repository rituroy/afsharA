dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
setwd(paste(dirSrc2,"AfsharA",sep=""))


##############################################

load("tmp.RData")

###########################################################
library(marray)
source(paste(dirSrc,"functions/heatmap.5.6.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.3.R",sep=""))

#source(paste("code/heatmap.5.6.R",sep=""))

#clId=rep(T,nrow(expr))

lineFlag=T

colList=c("skyblue","blue","yellow","purple","red")
colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
colList=c("blue","forestgreen","yellow","red")
colList2=c("skyblue","blue")
colList2=colList
colList2=c("blue","red")
colHM=c("red","blue","grey")
colHM=list(c("grey"),NULL,NULL)
colHM=c("grey","white","white")
colHM=list(c("grey","red","blue","orange","cyan","yellow","indianred","yellow2","skyblue","bisque3","indianred4"),NULL,NULL)
colHM=list(c("grey40","grey70","red","blue","orange","cyan","yellow","indianred","yellow2","skyblue","bisque3","indianred4"),NULL,NULL)

distMethod="pearson"
distMethod="cosine"
linkMethod="ward.D2"

outFormat="pdf"
outFormat=""
outFormat="png"

altTypeUniq1=1:2
altTypeUniq2=cbind(altTypeUniq1,as.character(1:length(altTypeUniq1)))

datadir=""

centrFlag="_noCentering"
centrFlag=""

subsetFlag=""

numPr=500
numPr=2000

pThres=10^-8
pThres=10^-6
pThres=0.05

compList=paste("_rnd",numPr,sep="")
compList=paste("_topVar",numPr,sep="")
#compList=paste("_",unique(ann$funcGrp),sep="")
compList=""

datFlag="_combatAdj"
datFlag=""

colGeneId="affyId"; colIdPV="FDR"; colNamePV="QV"

varList=c("monosomy3Type","Chrom1pLoss","Chrom6pgainBi","Chrom8qgain","GNAQ","GNA11","BAP1.loss","EIF1Ax","SF3B1")
datObj=list(expr=t(as.matrix(clin[,varList])),ann=data.frame(id=varList,geneId=varList,geneName=varInfo$varNameLong[match(varList,varInfo$varList)],stringsAsFactors=F),phen=clin[,c("id","CCGL.Pt.Number","Largestbasaldiam","Thickness","tnm4","path3","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","GEP3")])
for (k in c("Ciliarybodyinvolvment","EOE","Epithelioidcellsany")) {
    datObj$phen[,k]=as.character(datObj$phen[,k])
    datObj$phen[which(datObj$phen[,k]=="0"),k]="absent"
    datObj$phen[which(datObj$phen[,k]=="1"),k]="present"
}

tblCC=NULL
for (compFlag in compList) {
    compFName=sapply(compFlag,function(x) {y= strsplit(x," ")[[1]]; y=paste(y[1:min(c(length(y),3))],collapse="_"); gsub("-|/","_",y)}, USE.NAMES=F)
    if (length(grep("_rnd",compFlag))==1) {
        rndVec=paste("_rnd",1:4,sep="")
        #		rndVec=paste("_rnd",1:20,sep="")
    } else {
        rndVec=""
    }
    
    for (rndId in rndVec) {
        limFCmmu=c(-6,6)
        if (compFlag%in%c("_topSignif")) {
            load(paste(datadir,"dge_Homo_sapiens_",compFlag,".RData",sep=""))
            switch(compFlag,
            "_topSignif"={
                stat_1=stat1_4
            }
            )
            timeThis=as.integer(sub("hrs","",strsplit(compFlag,"_")[[1]][2]))
            compName1=paste(timeThis,"hr: TGFbeta vs. untreated",sep="")
        } else {
            if (substr(compFlag,1,nchar("_topVar"))=="_topVar") {
                compName1=paste("Top ",sub("_topVar","",compFlag)," most variable probesets",sep="")
            }
            if (substr(compFlag,1,nchar("_rnd"))=="_rnd") {
                compName1=paste("Random ",sub("_rnd","",compFlag)," probesets",sep="")
            } else {
                compName1=sub("_","",compFlag)
            }
        }
        for (transFlag in c("")) {
            if (F) {
                if (transFlag=="") {
                    subsetFlag=subsetName=""
                } else {
                    subsetFlag=paste("_",tolower(transFlag),sep="")
                    subsetName=paste(", ",transFlag,sep="")
                }
            }
            for (subsetFlag in c("")) {
                compName2=compName1
                if (subsetFlag=="") {
                    subsetName=""
                    prId=NULL
                    samId=1:nrow(datObj$phen)
                    sampleBar="cluster"
                    geneBar="clusterPr"
                    nClust=c(2,3)
                } else {
                    grpUniq=sub("_","",subsetFlag)
                    subsetName=paste(", ",grpUniq,sep="")
                    if (grpUniq=="gamma") grpUniq="_"
                    fName2=paste(datadir,"clusterInfoFeature",fName1,compFlag,centrFlag,rndId,".txt",sep="")
                    fName2=paste(datadir,"clusterInfoFeature",fName1,centrFlag,rndId,".txt",sep="")
                    prId=read.table(file=fName2, header=T, sep="\t", quote="", comment.char="", as.is=T)
                    prId=prId[,"probesetid"]
                    samId=which(tolower(datObj$phen$type3)==grpUniq)
                    sampleBar=""
                    geneBar=""
                    sampleBar="cluster"
                    geneBar="clusterPr"
                    nClust=c(NA,NA)
                }
                if (compFlag%in%c("_topSignif")) {
                    i1=which(stat_1[,colIdPV]<pThres)
                    if (length(i1)==0) next
                    fNameOut=paste(fName1,compFName,subsetFlag,centrFlag,"_",colNamePV,pThres,datFlag,sep="")
                    header=paste(compName2,subsetName,", ",colNamePV,"<",pThres,sep="")
                    dat0=eset$expr
                    switch(datFlag,
                    "_combatAdj"={dat0=exprCom
                    }
                    )
                    dat0=dat0[match(stat_1[,colGeneId][i1],rownames(dat0)),]
                } else {
                    #fNameOut=paste(fName1,compFName,subsetFlag,centrFlag,rndId,datFlag,sep="")
                    fNameOut=""
                    header=paste(compName2,subsetName,sep="")
                }
                expr=datObj$expr[,samId]
                annRow=datObj$ann
                phen=datObj$phen[samId,]
                phen$id2=sapply(phen$id,function(x) {strsplit(x,"_")[[1]][2]},USE.NAMES=F)
                
                annRowAll=annRow

                i2=1:nrow(expr)
                if (length(grep("_rnd",compFlag))==1) {
                    if (compFlag==paste("_rnd",numPr,sep="")) {
                        i1=1:numPr
                    }
                    header=paste(header,": ",length(i1)," random probesets",sep="")
                    geneBar="clusterPr"
                    #set.seed(5453)
                    i2=sample(1:nrow(expr),length(i1),replace=F)
                    expr=expr[i2,]
                    annRow=annRow[i2,]
                    i=1:nrow(expr)
                } else if (length(grep("_topVar",compFlag))==1) {
                    if (is.null(prId)) {
                        varGene=apply(expr,1,var,na.rm=T)
                        i2=order(varGene,decreasing=T)[1:numPr]
                    } else {
                        compName2=paste(compName2," based on all samples",sep="")
                        i2=match(prId,annRow[,"probesetid"])
                    }
                    header=paste(compName2,subsetName,sep="")
                    expr=expr[i2,]
                    annRow=annRow[i2,]
                    i=1:nrow(expr)
                } else if (length(grep("_top",compFlag))==1) {
                    header=paste(header,", n=",nrow(expr),sep="")
                    geneBar=""
                    geneBar="clusterPr"
                    i=1:nrow(expr)
                } else {
                    geneBar=""
                    geneBar="clusterPr"
                    expr=expr[i2,]
                    annRow=annRow[i2,]
                    i=1:nrow(expr)
                }
                if (!is.na(nClust[1])) {
                    nClust[1]=min(c(nrow(expr)-1,nClust[1]))
                    if (nrow(expr)<5) nClust[1]=NA
                }
                
                if (transFlag=="") {
                    j=1:ncol(expr)
                } else {
                    j=which(phen$translocation==transFlag)
                }
                
                
                arrayData=expr[i,j]
                annRow=annRow[i,]
                annCol=phen[j,]
                
                annColAll=datObj$phen
                annColAll$id2=annColAll$id
                
                if (F) {
                    if (centrFlag=="") {
                        centr=apply(arrayData,1,median,na.rm=T)
                        arrayData=arrayData-centr
                    }
                }
                
                varFList=varFName=NULL
                varFListAll=varFList
                varFNameAll=varFName
                
                varList=c("Largestbasaldiam","Thickness","tnm4","path3","Ciliarybodyinvolvment","EOE","Epithelioidcellsany","GEP3")
                varList=c("Largestbasaldiam","Thickness","tnm4","Ciliarybodyinvolvment","EOE","Epithelioidcellsany")
                #varName=paste(varList," ",sep="")
                #varName=paste(c("Largestbasaldiam","Thickness","TNM","Pathology","CiliaryBodyInvolvment","EOE","EpithelioidCell","GEP")," ",sep="")
                varName=paste(varInfo$varNameLong[match(varList,varInfo$varList)]," ",sep="")
                k=which(varList%in%names(annCol))
                varListAll=varList
                varNameAll=varName
                varList=varList[k]
                varName=varName[k]
                
                #nameRow=rep("",nrow(annRow))
                nameRow=annRow$geneName
                if (is.null(varFList)) {
                    colRow=NULL
                } else {
                    colRow=matrix(nrow=length(varFList),ncol=nrow(annRow))
                    for (varId in 1:length(varFList)) {
                        x=annRowAll[,varFList[varId]]
                        x[x==""]=NA; x=as.integer(as.factor(x))
                        grpUniq=sort(unique(x))
                        x=x[match(annRow$affyId,annRowAll$affyId)]
                        if (length(grpUniq)<=length(colList2)) {
                            colRow[varId,]=colList2[x]
                        } else if (length(grpUniq)<=length(colList)) {
                            colRow[varId,]=colList[x]
                        } else {
                            colRow[varId,]=rainbow(length(grpUniq))[x]
                        }
                    }
                    rownames(colRow)=varFName
                }
                
                nameCol=annCol$id2
                nameCol=rep("",nrow(annCol))
                colCol=NULL
                colCol=matrix(nrow=length(varList),ncol=nrow(annCol))
                for (varId in 1:length(varList)) {
                    if (varList[varId]%in%c("Largestbasaldiam","Thickness")) {
                        j=match(annCol$id,annColAll$id)
                        x=round(annColAll[,varList[varId]])
                        lim=range(x,na.rm=T)
                        #lim=quantile(x,probs=c(.1,.9),na.rm=T)
                        x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                        grpUniq=lim[1]:lim[2]
                        colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        colCol[varId,]=colColUniq[x[j]]
                    } else {
                        if (varList[varId]%in%c("time")) {
                            x=annColAll[,varList[varId]]
                        } else {
                            x=as.character(annColAll[,varList[varId]])
                        }
                        x[x==""]=NA; x=as.integer(as.factor(x))
                        grpUniq=sort(unique(x))
                        x=x[match(annCol$id,annColAll$id)]
                        if (length(grpUniq)<=length(colList2)) {
                            colCol[varId,]=colList2[x]
                        } else if (length(grpUniq)<=length(colList)) {
                            colCol[varId,]=colList[x]
                        } else {
                            colCol[varId,]=rainbow(length(grpUniq))[x]
                        }
                    }
                }
                rownames(colCol)=varName
                
                if (sampleBar=="cluster") {
                    switch(distMethod,
                    "cosine"={
                        dat=arrayData
                        if (any(apply(dat,2,sum,na.rm=T)==0)) dat=dat+1
                        distMat=getCosineDist(t(dat))
                    }
                    )
                    clustC=hclust(distMat, method=linkMethod)
                }
                if (geneBar=="clusterPr") {
                    switch(distMethod,
                    "cosine"={
                        dat=arrayData
                        if (any(apply(dat,2,sum,na.rm=T)==0)) dat=dat+1
                        distMat=getCosineDist(dat)
                    }
                    )
                    clustR=hclust(distMat, method=linkMethod)
                }
                
                print("summary(range(c(arrayData),na.rm=T))")
                print(summary(range(c(arrayData),na.rm=T)))
                arrayData2=arrayData
                arrayData2[arrayData==0]=NA
                if (T) {
                    for (j in 1:ncol(arrayData2)) {
                        for (k in 1:nrow(altTypeUniq2)) {
                            arrayData2[which(arrayData[,j]==as.integer(altTypeUniq2[k,2])),j]=k
                        }
                    }
                }
                limit1=c(0,1)
                if (is.null(colHM[[2]])) {
                    y=unique(c(arrayData2)); y=sort(y[!is.na(y)])
                    colHM[[1]]=colHM[[1]][1:max(y)]
                    if (length(colHM[[1]])==1) colHM[[1]]=rep(colHM[[1]],2)
                    limit1[2]=max(c(arrayData2),na.rm=T)
                }
                lim=150
                cexThis=c(2,2,0.5)
                main=NULL
                main=header
                cexThis=c(2,2,0.5)
                addTextThis=NULL
                if (lineFlag) {
                    lim=150
                    if (max(dim(arrayData))>lim) {
                        lineList=list(row=NULL,col=NULL,color=NULL)
                    } else {
                        lineList=list(row=seq(0,(nrow(arrayData2)))+0.5,col=seq(0,(ncol(arrayData2)))+0.5,color="grey")
                    }
                } else {
                    lineList=list(row=NULL,col=NULL,color=NULL)
                }
                
                subDir=""
                subDir=paste(sub("_","",compFName),"/",sep="")
                subDir="heatmap/"
                if (subDir!="" & !file.exists(subDir)){
                    dir.create(file.path(subDir))
                }
                if (outFormat=="png") {
                    margins=c(10,20)
                    margins=c(6,0.5)
                    margins=c(6,6)
                    png(paste(subDir,"heatmap",fNameOut,".png",sep=""),width=480*2,height=480*2)
                } else {
                    margins=c(12,5)
                    pdf(paste(subDir,"heatmap",fNameOut,".pdf",sep=""))
                }
                totalC=ncol(arrayData)
                #hcc=heatmap3(x=arrayData2, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=colRow, labCol=nameCol, labRow=nameRow, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit1, cexCol=cexThis[2], cexRow=cexThis[1], high=colHM[[1]], low=colHM[[2]], mid=colHM[[3]],lineRow=lineList$row, lineCol=lineList$col, lineColor=lineList$color, addText=addTextThis, cexText=cexThis[3])
                hcc=heatmap3(x=arrayData2, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=colRow, labCol=nameCol, labRow=nameRow, ncr=NA, ncc=NA, scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit1, cexCol=cexThis[2], cexRow=cexThis[1], high=colHM[[1]], low=colHM[[2]], mid=colHM[[3]],lineRow=lineList$row, lineCol=lineList$col, lineColor=lineList$color, addText=addTextThis, cexText=cexThis[3])
                dev.off()
                
                if (is.na(nClust[1])) {
                    clustId=paste("cluster",1,sep="")
                    i=1:nrow(annRow)
                } else {
                    clustId=cutree(clustR,k=nClust[1])[clustR$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    i=clustR$order
                }
                tbl=cbind(geneMut=annRow$geneName[i],clustId,order=1:nrow(annRow))
                write.table(tbl, paste(subDir,"clusterInfoFeature",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                
                if (is.na(nClust[2])) {
                    clustId=paste("cluster",1,sep="")
                    j=1:nrow(annCol)
                } else {
                    if (F) {
                        pdf(paste(subDir,"clusterSamples",fNameOut,".pdf",sep=""))
                        plot(clustC,main=paste("Sample clusters with ",nClust[2]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustC,k=nClust[2])
                        dev.off()
                    }
                    clustId=cutree(clustC,k=nClust[2])[clustC$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    j=clustC$order
                }
                tbl=cbind(annCol[j,which(!names(annCol)%in%c("id2","order"))],clustId,order=1:nrow(annCol))
                write.table(tbl, paste(subDir,"clusterInfoSample",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            }
        }
    }
}
subDir=""
subDir="legend/"
if (subDir!="" & !file.exists(subDir)){
    dir.create(file.path(subDir))
}
if (!is.null(colRow)) {
    for (varId in 1:length(varFListAll)) {
        if (length(grep("signif_",varFListAll[varId]))==1) {
            nm=""
            header=""
        } else {
            nm=varFListAll[varId]
            header=varFNameAll[varId]
        }
        if (outFormat=="png") {
            png(paste(subDir,"heatmapGeneColorBarLegend_",nm,".png",sep=""))
        } else {
            pdf(paste(subDir,"heatmapGeneColorBarLegend_",nm,".pdf",sep=""))
        }
        x=annRowAll[,varFListAll[varId]]
        grpUniq=table(x)
        grpUniq=names(grpUniq)
        k=1:length(grpUniq)
        if (length(grpUniq)<=length(colList2)) {
            sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=header)
        } else if (length(grpUniq)<=length(colList)) {
            sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=header)
        } else {
            sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=header)
        }
        dev.off()
    }
}
if (outFormat=="") {
    png("heatmapLegend.png",width=4*480,height=2*480)
    par(mfrow=c(2,4))
}
if (!is.null(colCol)) {
    for (varId in 1:length(varListAll)) {
        if (varListAll[varId]%in%c("Largestbasaldiam","Thickness")) {
            if (outFormat=="png") {
                png(paste(subDir,"heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""),width=480,height=0.3*480)
            } else if (outFormat=="pdf") {
                pdf(paste(subDir,"heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
            }
            x=round(annColAll[,varListAll[varId]])
            lim=range(x,na.rm=T)
            #lim=quantile(x,probs=c(.1,.9),na.rm=T)
            grpUniq=lim[1]:lim[2]
            colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
            heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1)],median(1:length(colColUniq))),main=varNameAll[varId])
        } else {
            if (outFormat=="png") {
                png(paste(subDir,"heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""),width=480,height=480)
            } else if (outFormat=="pdf") {
                pdf(paste(subDir,"heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
            }
            if (varList[varId]%in%c("time")) {
                x=annColAll[,varListAll[varId]]
            } else {
                x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
            }
            grpUniq=table(x)
            #grpUniq=paste(names(grpUniq)," (",grpUniq,")",sep="")
            grpUniq=names(grpUniq)
            k=1:length(grpUniq)
            if (length(grpUniq)<=length(colList2)) {
                sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=varNameAll[varId],cex=1.5)
            } else if (length(grpUniq)<=length(colList)) {
                sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=varNameAll[varId],cex=1.5)
            } else {
                sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=varNameAll[varId],cex=1.5)
            }
        }
        if (outFormat!="") dev.off()
    }
}
colHMAll=colHM
fNameOut2=""
if (length(colHMAll[[1]])==1) {
    if (outFormat=="png") {
        png(paste(subDir,"heatmapColorRange",fNameOut2,".png",sep=""),width=480,height=140)
    } else if (outFormat=="pdf") {
        pdf(paste(subDir,"heatmapColorRange",fNameOut2,".pdf",sep=""))
    }
    if (F) {
        heatmapColorBar=function(limit,cols=c("green","red","black")) {
            try <- maPalette(high=cols[1], low=cols[2], mid=cols[3],k=15)
            maColorBar(try, scale=limit,k=3)
        }
    }
    heatmapColorBar(limit=limit1,cols=unlist(colHMAll),main="Alteration")
    dev.off()
} else {
    width = 480; height = 480
    if (outFormat=="png") {
        png(paste(subDir,"heatmapColorLegend",fNameOut2,".png",sep=""),width=width,height=height)
    } else if (outFormat=="pdf") {
        pdf(paste(subDir,"heatmapColorLegend",fNameOut2,".pdf",sep=""))
    }
    grpUniq=sort(altTypeUniq2[,1])
    grpUniq=c("mutation / monosomy 3","partial monosomy 3")
    cexThis=NULL
    cexThis=1.5
    if (outFormat=="pdf") cexThis=1
    sampleColorLegend(tls=grpUniq,col=colHMAll[[1]],legendTitle="Mutation status",cex=cexThis)
    if (outFormat!="") dev.off()
    
}
if (outFormat=="") dev.off()

###########################################################
###########################################################