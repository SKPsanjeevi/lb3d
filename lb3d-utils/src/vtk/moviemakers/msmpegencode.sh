mencoder mf://\*.jpg -mf w=384:h=288:fps=25:type=jpg -ovc lavc -lavcopts vcodec=msmpeg4:mbd=2:vbitrate=3000:keyint=250 -o msoutput.avi
