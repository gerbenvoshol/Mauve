#
# read a list of bp localization values and make plots
#

rundesc <- read.table( "run_description.txt" )
a <- read.table("scorelist.txt");
b <- a + 1
c <- log(b)
z<-as.matrix(c)

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


postscript( "scoreplot.ps",onefile=TRUE,paper="special",bg="transparent",width=8,height=8,horizontal=FALSE )
image(seq(0,rundesc[1,3]*(rundesc[1,1]-1),,ncol(z)),seq(0,rundesc[1,4]*(rundesc[1,2]-1),,nrow(z)),t(z),,,, rev(mycolors))

legmat <- seq(1,4096,,4096);
legmat <- as.matrix(legmat);
postscript( "bp_legend.ps",onefile=TRUE,paper="special",bg="transparent",width=5.5,height=2,horizontal=FALSE )
image(seq(min(z),max(z),,4096),seq(1,1,,1),legmat,,,, rev(mycolors), xlab="log10 of absolute distance between predicted and correct bp")

