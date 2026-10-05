# Display the positive-intensity counts from All_expression.xlsx.
# Run from proteomics/. Counts are supplied in Protein_detection_counts.csv;
# the main QC script also derives these counts directly from the workbook.
out <- "Supplementary_Figures"
dir.create(out,showWarnings=FALSE)
d <- read.csv("Protein_detection_counts.csv")
fences <- sapply(c("F1","F2"),function(g){x<-d$Proteins_Identified[d$Group==g];unname(quantile(x,.25)-1.5*IQR(x))})
stopifnot(identical(as.numeric(round(fences)),c(2377,1446)))
draw <- function(){
 par(family="Arial",oma=c(1.5,3,4.8,.5),mfrow=c(2,1),mar=c(2.3,2.8,2.1,.6),mgp=c(2.2,.65,0),las=1)
 for(g in c("F1","F2")){
  s<-d[d$Group==g,];s<-s[order(s$Proteins_Identified),]
  plot(NA,xlim=c(0,4000),ylim=c(.4,6.6),xaxs="i",yaxs="i",axes=FALSE,xlab="",ylab="")
  color<-if(g=="F1")"#EA6E5B" else "#6CC6D8"
  rect(0,seq_len(nrow(s))-.43,s$Proteins_Identified,seq_len(nrow(s))+.43,col=color,border="black",lwd=1.3)
  axis(2,at=seq_len(nrow(s)),labels=s$Mouse_ID,tcl=-.25,cex.axis=1.0)
  if(g=="F2")axis(1,at=seq(0,4000,1000),tcl=-.25,cex.axis=1.0)
  segments(0,.4,0,6.6,lwd=1.3)
  if(g=="F2")segments(0,.4,4000,.4,lwd=1.1)
  text(s$Proteins_Identified+26,seq_len(nrow(s)),s$Proteins_Identified,adj=0,cex=.88)
  abline(v=fences[g],lty=2,lwd=1.5)
  text(fences[g]-50,3.5,paste("Lower IQR fence:",round(fences[g])),srt=90,cex=.93,font=3)
  par(xpd=NA)
  rect(0,6.68,4000,7.24,col="#E6E6E6",border="#4A4A4A",lwd=.9)
  text(2000,6.96,g,font=2,cex=1.02)
  par(xpd=FALSE)
 }
 mtext("QC: Proteomic Detection & Outlier Analysis",side=3,outer=TRUE,line=3,adj=.05,cex=1.4)
 mtext("Dashed lines mark the lower IQR fence (Q1 - 1.5 x IQR) for each generation.",side=3,outer=TRUE,line=1.5,adj=.05,cex=1.03)
 mtext("Mouse ID",side=2,outer=TRUE,line=1.35,las=0,cex=1.02)
 mtext("Number of Proteins Detected (>0)",side=1,outer=TRUE,line=.1,cex=1.03)
}
tiff(file.path(out,"SuppFig_QC_ProteinCounts_Stats.tiff"),width=3600,height=2100,res=300,compression="lzw",type="cairo");draw();dev.off()
