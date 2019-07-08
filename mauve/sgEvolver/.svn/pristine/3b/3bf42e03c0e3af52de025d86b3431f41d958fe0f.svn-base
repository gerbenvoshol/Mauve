# reads in data from a file stored at scorelist.txt
# format: m x n matrix of scores
# reads a run description from a file called run_description.txt
# format: num_cols num_rows x_increment y_increment
#       : x axis name
#       : y axis name
#       : title
# writes out a graph to a file called scoreplot.ps
#
# (c) 2003-2005 Aaron Darling

# read and format the data
mscores <- read.table( "scorelist.txt" )
rundesc <- read.table( "run_description.txt" )
matvals <- matrix( mscores[,1], rundesc[1,1], rundesc[1,2] )

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
postscript( "scoreplot.ps",onefile=TRUE,paper="special",bg="transparent",width=8,height=8,horizontal=FALSE )
image(seq(0,,rundesc[1,3],rundesc[1,1]),seq(0,,rundesc[1,4],rundesc[1,2]),matvals,c(0,1),,, mycolors )
title( "Alignment accuracy plot" )
dev.off()
