## Oldrich Rakovec 
## 8 March 2013
## contact (from 1dec 2013: oldrich.rakovec ET ufz.de)
##
## This R scripts plots most important parameters of the FUSE-016 model structure of the
## work published by Rakovec et al. (2014, WRR, doi:10.1002/2013WR014063) colored by differerent
## parameter values;
## 
## ref: Rakovec, O., M. C. Hill, M. P. Clark, A. H. Weerts, A. J. Teuling, R. Uijlenhoet (2014), Distributed
## Evaluation of Local Sensitivity Analysis (DELSA), with application to hydrologic models,
## Water Resour. Res., 50, 1-18, doi:10.1002/2013WR014063.
## 
## #################################################################################################
## #################################################################################################

rm(list = ls(all = TRUE))
graphics.off()

## packages to be loaded
library(RColorBrewer)
library(ncdf4)
ncfunget=ncvar_get

## margins for all 4 parameters
margi2=matrix(NA,ncol=4, nrow=4) 
margi2[,1]=rep(seq(from=0.20, to=0.8, len=4))
margi2[,2]=rep(seq(from=0.20, to=0.8, len=4)+0.17)
margi2[,3]=c(rep(0.25,4))
margi2[,4]=c(rep(0.9,4))
print(margi2)

## selection of 1000 items for plotting purposes
NROW=1000
delsa_index= read.table("../data_input/fuse016_delsa_raw_rmse.csv",skip=0,header=TRUE,sep=",",nrow=NROW )
perf_index = unlist(read.table("../data_input/fuse016_raw_rmse.csv",skip=1 , nrow=NROW))
parval     = read.table("../data_input/fuse016_parameters_base.csv",skip=0,header=TRUE,sep=",",nrow=NROW )

PARNMS=c("MAXWATR_1","MAXWATR_2","FRACTEN", "PERCRTE","PERCEXP","AXV_BEXP","TIMEDELAY")

plotName=paste("figure_03.pdf")
pdf(plotName, width = 7, height = 2., family="Helvetica")
par(oma=c(0.,0.,0,0), mar=c(0.2,0.2,0.2,0.2), las=1)#, lwd=2, cex.lab=2.7, cex.axis=2.7)
par(mgp = c(1.35, 0.2, 0))
par(plt=c(0.12,0.88,0.22,0.9))
par(cex.main=0.8)
colo=rev(c(1:6,8))

selected=which(names(delsa_index)%in%PARNMS)
npars=length(selected)

ii = 0

parplot=c(1,3,6,7)

## color each param rmse combination given parameter valus
for (co in selected[parplot]){

    ii=ii+1
    plot.new()
    prvalco=parval[,co]
    
    zq=c(prvalco) ## addition for fixing scales
    leng=10 # number of colors in the color bar
    ## select colorful scheme:
    col_tab=rev(brewer.pal(leng,"Spectral"))
    Col1 <- col_tab[as.numeric(cut(zq,breaks = leng))]

    ## this is for two scatterplots with different extent but having the same color scale
    if(co == 3){
        Col1 = 1
    }
    
    mx=1
    
    for (pn1 in selected[parplot]){
        
        par(plt=margi2[mx,], new=TRUE)
        mx =mx+1
        delsa=delsa_index[,pn1]

        ## #############
        ## plot parameter value vs. model output:
        plot(delsa,perf_index, col=Col1, pch=20, cex=0.75,
             ylab=ifelse(mx==2,"raw_rmse",""), xlab=ifelse(mx==3,expression(italic(S[i]^L)),"") ,
             xaxt="n", yaxt="n", xlim=c(0,1),main="")
        title(names(delsa_index)[selected[parplot]][mx-1], font.main = 1, line = 0.15)
        axis(1, labels=c("0","0.5","1"), at=c(0,0.5,1),  tck=0.05)        
        axis(2, labels=ifelse(mx==2,TRUE,FALSE), tck=0.05)        

    } 

    if(co == 3){

        print("param 3")

    }else{

        ## plot the corresponding color scheme:
        lower=c(50,25,0.05,0.01,1,0.001,0.1)[parplot]
        upper=c(500,250,0.95,1000,20,3,2)[parplot]
        
        par(plt=c(0.08,0.11,0.20,0.9), new=TRUE)
        zz=array(1:leng,dim=c(1,leng))
        image(1,1:leng,zz, col = col_tab, yaxt="n",xlab="",  ylab="", xaxt = "n", cex.axis = 1.2, cex.lab = 1.5)
        labe=round(seq(lower[ii],upper[ii],len=leng+1),digits=1)
        axis(side=2, at= c(1:(leng+1))-0.5, labels=labe, tck=-0.05)
        box()

        par(plt=c(0.0,0.15,0.06,0.14), new=TRUE)
        plot(-999,-999, xlab="", ylab="", xaxt="n", yaxt="n", axes=FALSE, ylim=c(0,2), xlim=c(0,2))
        text(1,1, names(delsa_index)[co] , cex=1)

    }

}
graphics.off()
