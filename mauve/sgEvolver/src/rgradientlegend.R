# Generates a legend for sgEvolver accuracy plots
# Writes the legend to the file "legendplot.ps" in the current directory
# process this file using `R CMD BATCH rgradientlegent.R`
#
# (c) 2003-2005 Aaron Darling

#generate one row per color
matvals <- matrix( seq( 1, 3072 ) / 3072, 1, 3072 );

# generate a color table
# first do black to red
hues <- seq( 0, 0,, 1024 );
sats <- seq( 1, 1,, 1024 );
values <- seq( 0, 1,, 1024 );
alphas <- seq( 1, 1,, 1024 );
# then do red to yellow
hues <- c(hues,seq( 0, (1/6),, 1024 ));
sats <- c(sats,seq( 1, 1,, 1024 ));
values <- c(values,seq( 1, 1,, 1024 ));
alphas <- c(alphas,seq( 1, 1,, 1024 ));
# then do yellow to white
hues <- c(hues,seq( (1/6), (1/6),, 1024 ));
sats <- c(sats,seq( 1, 0,, 1024 ));
values <- c(values,seq( 1, 1,, 1024 ));
alphas <- c(alphas,seq( 1, 1,, 1024 ));
mycolors = hsv( hues, sats, values, alphas );

# write to a postscript file
postscript( "legendplot.ps",onefile=TRUE,paper="special",bg="transparent",width=2,height=5.5,horizontal=FALSE )
par(xaxt="n")
image( z=matvals, col=mycolors )
#zlim=c(0,1)
#image(seq(0,,rundesc[1,3],rundesc[1,1]),seq(0,,rundesc[1,4],rundesc[1,2]),matvals,c(0,1),c(0,rundesc[1,1]*rundesc[1,3]),c(0,rundesc[1,2]*rundesc[1,4]), mycolors )
#"Legend"
title( ,,, "Alignment Accuracy")
dev.off()
