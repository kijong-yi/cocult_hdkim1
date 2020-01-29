
cat << EOF 1>&2
--------------
migec pipeline
--------------
v1.1 seperately checkout A, B
`date`

args: $@

EOF

. ~kjyi/src/parse2 $@ << EOF
-o		odir	"."	output directory
--a1	a1		""	master barcode sequence for TCRA
--a2	a2		""	slave barcode sequence for TCRA
--b1	b1		""	master barcode sequence for TCRB
--b2	b2		""	slave barcode sequence for TCRB
--f1	f1		""	input fastq 1
--f2	f2		""	input fastq 2
EOF


if [ -f $f1 ] && [ -f $f2 ] ; then
	echo Input files present: $f1, $f2
else
	echo Input files are not exist;	exit 1
fi

mkdir -p $odir
cat << EOF > $odir/barcodes_A.txt
TCRA	$a1	$a2
EOF
cat << EOF > $odir/barcodes_B.txt
TCRB	$b1	$b2
EOF

echo "[[ Checkout ]]"
mkdir -p $odir/checkout/A
mkdir -p $odir/checkout/B
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Checkout -cute $odir/barcodes_A.txt $f1 $f2 $odir/checkout/A
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Checkout -cute $odir/barcodes_B.txt $f1 $f2 $odir/checkout/B

echo "[[ Histogram ]]"
mkdir -p $odir/histogram/A
mkdir -p $odir/histogram/B
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Histogram $odir/checkout/A $odir/histogram/A
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Histogram $odir/checkout/B $odir/histogram/B

pushd $odir/histogram/A
R --slave << EOF

#
# Plots a fancy oversequencing histogram from tables pre-computed by Histogram util
#
# usage:
# \$cd histogram_output_dir/
# \$RScript histogram.R
#

require(ggplot2); require(reshape)

build_df <- function(stat, units, sweep_df = NULL) {
	df<- read.table(paste(stat, units, ".txt", sep= ""), header=FALSE, comment="", sep = "\\t")
  x<-df[1,3:ncol(df)]
	id<-df[2:nrow(df), 1] 
	if (is.null(sweep_df)) {
		sweep_df <- df
	}
	df[2:nrow(df),3:ncol(df)] <- sweep(df[2:nrow(df),3:ncol(df)], 1, rowSums(sweep_df[2:nrow(df),3:ncol(df)]), FUN = "/")
	      
	df.m <- melt(df[2:nrow(df), c(1,3:ncol(df))])
	df.m\$x <- as.vector(t(x[df.m\$variable]))
	df.m\$s <- rep(stat, nrow(df.m))
	list(sweep_df, df.m)						   
}
   
for (units in c("", "-units")) {  
   dfo <- build_df("overseq", units)
   dfc <- build_df("collision1", units, dfo[[1]])     
   df <- rbind(dfo[[2]], dfc[[2]])
   
   pdf(paste("histogram", units, ".pdf", sep= ""))
   
	 print(
		ggplot(df, aes(x=x,y=value))+geom_smooth()+
			scale_x_log10(name = "MIG size, reads", expand=c(0,0), limits=c(1, 10000), breaks = c(1:10,100,1000,10000), labels = c("1", rep("", 8), 10, 100, 1000, 10000), oob=scales::rescale_none)+
			scale_y_continuous(name = "",expand=c(0,0), limits = c(0, max(df\$value)), oob=scales::rescale_none)+theme_bw()+theme(panel.grid.minor = element_blank()) + facet_grid(s~.)
	 )
   
   dev.off()  
}
EOF

R --slave << EOF

#
# Plots position-weight matrices for UMI sequences based on data precomputed by Histogram util
#
# usage:
# \$cd histogram_output_dir/
# \$RScript pwm.R
#

require(seqLogo) # available @ bioconductor

# read in
logo <- function(prefix) {
	df<-read.table(paste(prefix,".txt",sep=""))
	rownames(df)<-df[,1]
	df<-df[,2:ncol(df)]
	df[, ] <- apply(df[, ], 2, as.numeric)
	
	# build seqlogo
	df.p <- makePWM(df)
	pdf(paste(prefix,".pdf",sep=""))
	seqLogo(df.p, ic.scale = F)
	dev.off()
	
	con <- file(paste(prefix,"-stats.txt",sep=""))
	sink(con, append=TRUE)
	sink(con, append=TRUE, type="output")
	
	# stats
	print(paste("IC =",sum(df.p@ic)))
	
	#entropy
	h<-sum(2-df.p@ic)
	print(paste("H =",h))
	
	#correlation with pos
	print(paste("R(Hi,i) =", cor(-df.p@ic,1:length(df.p@ic))))
	
	#number of variants
	print(paste("N_obs =",2^h))
	
	#theoretical number of variants
	print(paste("N_exp =",2^(2*ncol(df))))
	sink()
}

logo("pwm-summary")
logo("pwm-summary-units")
EOF

popd

pushd $odir/histogram/B
R --slave << EOF

#
# Plots a fancy oversequencing histogram from tables pre-computed by Histogram util
#
# usage:
# \$cd histogram_output_dir/
# \$RScript histogram.R
#

require(ggplot2); require(reshape)

build_df <- function(stat, units, sweep_df = NULL) {
	df<- read.table(paste(stat, units, ".txt", sep= ""), header=FALSE, comment="", sep = "\\t")
  x<-df[1,3:ncol(df)]
	id<-df[2:nrow(df), 1] 
	if (is.null(sweep_df)) {
		sweep_df <- df
	}
	df[2:nrow(df),3:ncol(df)] <- sweep(df[2:nrow(df),3:ncol(df)], 1, rowSums(sweep_df[2:nrow(df),3:ncol(df)]), FUN = "/")
	      
	df.m <- melt(df[2:nrow(df), c(1,3:ncol(df))])
	df.m\$x <- as.vector(t(x[df.m\$variable]))
	df.m\$s <- rep(stat, nrow(df.m))
	list(sweep_df, df.m)						   
}
   
for (units in c("", "-units")) {  
   dfo <- build_df("overseq", units)
   dfc <- build_df("collision1", units, dfo[[1]])     
   df <- rbind(dfo[[2]], dfc[[2]])
   
   pdf(paste("histogram", units, ".pdf", sep= ""))
   
	 print(
		ggplot(df, aes(x=x,y=value))+geom_smooth()+
			scale_x_log10(name = "MIG size, reads", expand=c(0,0), limits=c(1, 10000), breaks = c(1:10,100,1000,10000), labels = c("1", rep("", 8), 10, 100, 1000, 10000), oob=scales::rescale_none)+
			scale_y_continuous(name = "",expand=c(0,0), limits = c(0, max(df\$value)), oob=scales::rescale_none)+theme_bw()+theme(panel.grid.minor = element_blank()) + facet_grid(s~.)
	 )
   
   dev.off()  
}
EOF

R --slave << EOF

#
# Plots position-weight matrices for UMI sequences based on data precomputed by Histogram util
#
# usage:
# \$cd histogram_output_dir/
# \$RScript pwm.R
#

require(seqLogo) # available @ bioconductor

# read in
logo <- function(prefix) {
	df<-read.table(paste(prefix,".txt",sep=""))
	rownames(df)<-df[,1]
	df<-df[,2:ncol(df)]
	df[, ] <- apply(df[, ], 2, as.numeric)
	
	# build seqlogo
	df.p <- makePWM(df)
	pdf(paste(prefix,".pdf",sep=""))
	seqLogo(df.p, ic.scale = F)
	dev.off()
	
	con <- file(paste(prefix,"-stats.txt",sep=""))
	sink(con, append=TRUE)
	sink(con, append=TRUE, type="output")
	
	# stats
	print(paste("IC =",sum(df.p@ic)))
	
	#entropy
	h<-sum(2-df.p@ic)
	print(paste("H =",h))
	
	#correlation with pos
	print(paste("R(Hi,i) =", cor(-df.p@ic,1:length(df.p@ic))))
	
	#number of variants
	print(paste("N_obs =",2^h))
	
	#theoretical number of variants
	print(paste("N_exp =",2^(2*ncol(df))))
	sink()
}

logo("pwm-summary")
logo("pwm-summary-units")
EOF

popd





echo "[[ Assemble ]]"
mkdir -p $odir/assemble/A
mkdir -p $odir/assemble/B
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar AssembleBatch --force-collision-filter --force-overseq 2 $odir/checkout/A $odir/histogram/A $odir/assemble/A
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar AssembleBatch --force-collision-filter --force-overseq 2 $odir/checkout/B $odir/histogram/B $odir/assemble/B

echo "[[ CdrBlast ]]"
mkdir -p $odir/cdrblast/A
mkdir -p $odir/cdrblast/B
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar CdrBlastBatch -R TRA $odir/checkout/A $odir/assemble/A $odir/cdrblast/A
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar CdrBlastBatch -R TRB $odir/checkout/B $odir/assemble/B $odir/cdrblast/B

echo "[[ filter ]]"
mkdir -p $odir/cdr_filter/A
mkdir -p $odir/cdr_filter/B
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar FilterCdrBlastResultsBatch $odir/cdrblast/A $odir/cdr_filter/A
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar FilterCdrBlastResultsBatch $odir/cdrblast/B $odir/cdr_filter/B


echo "[[ reporting ]]"
mkdir -p $odir/report/A
mkdir -p $odir/report/B
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Report -c $odir/checkout/A -i $odir/histogram/A -a $odir/assemble/A -b $odir/cdrblast/A -f $odir/cdr_filter/A $odir/report/A
java -jar /home/users/kjyi/tools/migec/migec-1.2.9.jar Report -c $odir/checkout/B -i $odir/histogram/B -a $odir/assemble/B -b $odir/cdrblast/B -f $odir/cdr_filter/B $odir/report/B
echo "[[ editing report RMD script ]]"
sed -i '1{s,^,Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc");,}' $odir/report/A/*.R
sed -i '1{s,^,Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc");,}' $odir/report/B/*.R
sed -i '2{s!)!, knit_root_dir = "'$(pwd)'")!}' $odir/report/A/*.R
sed -i '2{s!)!, knit_root_dir = "'$(pwd)'")!}' $odir/report/B/*.R
sed -i '200,230{s,df$READS_DROPPED_WITHIN_MIG / df$READS_TOTAL,(df$READS_DROPPED_WITHIN_MIG_1 + df$READS_DROPPED_WITHIN_MIG_2) / df$READS_TOTAL,}' $odir/report/A/*.Rmd
sed -i '200,230{s,df$READS_DROPPED_WITHIN_MIG / df$READS_TOTAL,(df$READS_DROPPED_WITHIN_MIG_1 + df$READS_DROPPED_WITHIN_MIG_2) / df$READS_TOTAL,}' $odir/report/B/*.Rmd
sed -i '168{s!levels=\(.*\)!levels=unique(\1)!}' $odir/report/A/*.Rmd
sed -i '168{s!levels=\(.*\)!levels=unique(\1)!}' $odir/report/B/*.Rmd
echo "[[ running edited report RMD script ]]"
Rscript $odir/report/A/*.R
Rscript $odir/report/B/*.R

echo "[[ done ]]"



