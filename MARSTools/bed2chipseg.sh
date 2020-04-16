#!/usr/bin/env bash


function get_fasta(){

    #check the size of the input sequence of reduce the processing time form large files (assumes sorted file)

    lenbed=$(wc -l $bed_in | cut -f1 -d " ")
    lenuse=$(($lenbed/10))
    if [ "$lenbed" -gt 5000 ];
        then
            cut -f 1,2,3 $bed_in | head -$lenuse | bed-widen -width $len  >$bed_wide.bed
     else
            cut -f 1,2,3 $bed_in | bed-widen -width 100  >$bed_wide.bed
    fi

    #Extract negative bed and the fasta
    python extractnegative.py $bed_wide.bed $bed_wide.negbed 500

	fastaFromBed -tab -fi $hg -bed $bed_wide.bed -fo $bed_wide.fas
	fastaFromBed -tab -fi $hg -bed $bed_wide.negbed -fo $bed_wide.negfas

    cut -f 4 $bed_in >/tmp/f1
	cut -f 1 $bed_wide.fas >/tmp/f2
	cut -f 2 $bed_wide.fas >/tmp/f3
	paste /tmp/f2 /tmp/f1 /tmp/f3  >$bed_wide.fas

	cut -f 4 $bed_in >/tmp/f1
	cut -f 1 $bed_wide.negfas >/tmp/f2
	cut -f 2 $bed_wide.negfas >/tmp/f3
	paste /tmp/f2 /tmp/f1 /tmp/f3  >$bed_wide.negfas

	python removemasked.py $bed_wide.fas $bed_wide.fa
	python removemasked.py $bed_wide.negfas $bed_wide.negfa

    # use the length of the available sequences to determine the size of test and negative sequences

    lenfa=$(wc -l $bed_wide.fa | cut -f1 -d " ")
    len_negfa=$(wc -l $bed_wide.negfa | cut -f1 -d " ")
    lenuse=$(($lenbed/20))

    if [ $lenuse -gt 500 ];
        then
            lenuse=$lenuse
     else
            lenuse=500
     fi


    if [ "$lenfa" -gt $lenuse ] &&  [ "$len_negfa" -gt $lenuse ];
        then
            cutoff=$lenuse
    elif [ "$lenfa" -lt "$len_negfa" ] ;
        then
            cutoff=$lenfa
     else
            cutoff=$len_negfa
    fi

    head -$cutoff $bed_wide.fa >$bed_wide.posneg
	head -$cutoff $bed_wide.negfa >>$bed_wide.posneg

    #clean up the temporary files
	#rm $bed_wide.neg*
	#rm $bed_wide.fa*
	rm /tmp/f*
}
bed_in=$1
bed_wide=$2
hg=$3
len=${4:-100} #Use a default f 100 if not provided
get_fasta