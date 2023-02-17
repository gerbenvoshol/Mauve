# Mauve
Mauve is a system for constructing multiple genome alignments in the presence of large-scale evolutionary events such as rearrangement and inversion. Multiple genome alignments provide a basis for research into comparative genomics and the study of genome-wide evolutionary dynamics.

Mauve has been developed with the idea that a multiple genome aligner should require only modest computational resources. It employs algorithmic techniques that scale well in the lengths of sequences being aligned. For example, a pair of Y. pestis genomes can be aligned in under a minute, while a group of 9 divergent Enterobacterial genomes can be aligned in a few hours. However, the current algorithm’s compute time (progressiveMauve) scales cubically in the number of genomes to align, making it unsuitable for datasets containing more than 50-100 bacterial genomes.

Mauve development began at the University of Wisconsin-Madison with a team including Aaron Darling, Bob Mau, and Nicole Perna. Several others have contributed development to aspects of the Mauve software in the time since. Unfortunately, the software seems to be unmaintained at this time. 

This version has a patched mauve-gui which allows it to run on modern systems.

The orignal unpatched Mauve build instructions are included in the file `BUILD_INSTRUCTIONS.txt`. 

On `Ubuntu 19.04` it is possible to install Mauve using:

`sudo apt-get install mauve`

Otherwise download the repository and unzip it.

```
cd mauve/mauve
ant dist
```

the build is in the dist folder.

NOTE: for people with a high resolution display it is possible to scale the program by changing:

```
JAVA_ARGS="-Xms200M -Xmx500M"

into

JAVA_ARGS="-Xms200M -Xmx500M -Dsun.java2d.uiScale=2.5"
```
