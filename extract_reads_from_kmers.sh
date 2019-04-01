### reformat the kmer table for use in script kfastqfilter
### IN = the kmer table 

IN=$1
OUT=${IN}.ed
echo "# filtered kmer table" > $OUT
awk ' $2 >= 0 && $1 !~ /[#+N]/ ' $IN >> $OUT
COUNT=`grep -v "[#N]" $OUT|wc -l`
SUM=`awk '{ sum+=$2} END {print sum}' $OUT`
echo -e "+\t$SUM\t$COUNT\t31" >> $OUT

### extracting the reads using script kfastqfilter

kfastqfilter ${OUT} \
read1.fq read2.fq /dev/null /dev/null \
extracted_read1.fq extracted_read2.fq 1 

### the last number = how many kmers of the list the read must have at least
### see further explanations by executing kfastqfilter
