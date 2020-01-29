cd /home/users/kjyi/Projects/cocult/hdkim1

grep -v "^#" sample3.txt | while read num name am as bm bs r1 r2 ; do
	#cat << EOF | qsub -q day -l nodes=1:ppn=1 -e /dev/null -o /dev/null	
	cat << EOF > cmd/mm.$num.$name.sh	
cd /home/users/kjyi/Projects/cocult/hdkim1
mkdir -p migec/$name
sh /home/users/kjyi/Projects/cocult/hdkim1/pipe.sh \\
	-o migec/$name \\
	--a1 $am \\
	--a2 $as \\
	--b1 $bm \\
	--b2 $bs \\
	--f1 $r1 \\
	--f2 $r2 &> migec/$name/migec_pipeline.log
EOF
done
