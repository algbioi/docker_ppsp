#Author: Kaustubh R. Patil patil AT mpi-inf.mpg.de (c) 2010
#Thanks to Lars Steinbrueck for some tweaks.
#THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN  NO EVENT SHALL CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH  DAMAGE.

library("RSvgDevice")

scaffold_consistency_barplot <- function(file_base,min_l_scaf=0,main="",sub="",put_legend=TRUE,file_svg=NA) {

#options
#file_base: the base file name without _graph* (usually the predictions file name)
#min_l_scaf: what scaffolds to plot
#main and sub: The plot titles
#put_legend: Whether to put a legend
#file_svg: SVG output file (will be plotted on screen if this is NA)

mat_file <- paste(file_base,"_graph_mat.txt",sep="")
col_file <- paste(file_base,"_graph_col.txt",sep="")
rname_file <- paste(file_base,"_graph_rname.txt",sep="")
legend_file <- paste(file_base,"_graph_legend.txt",sep="")

t <- read.table(mat_file)
c <- as.vector(read.table(col_file)[,1])
r <- as.vector(as.matrix(read.table(rname_file)[,1]))
legend <- read.table(legend_file,sep="\t")
l <- as.vector(legend[,1])
f <- colors()[as.vector(legend[,2])]
m <- coordinate2matrix(t,0)

m <- t(m)

#NOTE!!! please check the following as the scaffold names (i.e. variable r) gets empty and hence no scaffold names printed
#remove the empty rows, this can happen if the files were generated using minimum scaffold length limit
csum <- apply(m,2,sum)
ii <- which(csum==0,arr.ind=TRUE)
if(length(ii)>0) {
	m <- m[,-ii]
	#r <- r[-ii]
}

if(min_l_scaf!=0) {
	csum <- apply(m,2,sum)
	ii <- which(csum<min_l_scaf,arr.ind=TRUE)
	if(length(ii)>0) {
        	m <- m[,-ii]
	        r <- r[-ii]
	}
}

d.order <- order(colSums(m),decreasing=F)
colnames(m) <- r

if(!is.na(file_svg)) 
	devSVG(file=file_svg)

#par(bg = "lightgrey")
par(pty="s")
p <- barplot(m[,d.order],col=colors()[c],horiz=T,main=main,sub=sub,xlab="Length",ylab="Scaffold identifier",las=1,border="white",cex.names=1)

if(put_legend)
	legend("bottomright",l,fill=f,inset=c(0,0.04),bty="n")

if(!is.na(file_svg))
	dev.off()

#ggplot(as.data.frame(m[,d.order])) + geom_bar()

}

#helper function
coordinate2matrix <- function(coord,clear) {
n_lines <- dim(coord)[1]

rows <- max(sapply(coord[,1], as.integer))
cols <- max(sapply(coord[,2], as.integer))
mat <- matrix(rep(clear,rows*cols),rows,cols)
for(i in 1:n_lines) {
	nrow <- as.integer(coord[i,1])
	ncol <- as.integer(coord[i,2])
	val <- as.double(coord[i,3])
	mat[nrow,ncol] <- val
}
mat
}
