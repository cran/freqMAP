`add.p.value.legend` <-
function(alpha1,alpha2,cex=1){
  
  legend(x="bottomleft",
         #x=xlims[1],y=ylims[1]+(ylims[2]-ylims[1])/lp,
         legend=paste("p <",c(alpha1,alpha2)),density=c(-1,-1),#density=c(28,56),
         fill=c("gray90","darkgray"),cex=cex,angle=c(90,90),bty="n")
  
}

