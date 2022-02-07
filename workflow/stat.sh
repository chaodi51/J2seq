
## stats of the mapping and dedup
rm -f mapped_reads.txt
touch -f mapped_reads.txt
for i in `sed 1d sample_table.tsv | cut -f1`;do
    echo $i; echo $i>name;
	samtools flagstat merged_bam/$i\_merged.bam > foo
	cat foo |grep "total"|head -1|cut -d " " -f1 > total
	cat foo |grep "mapped"|head -1|cut -d " " -f1 > uniq
	samtools flagstat merged_bam/$i\_merged_dedup_sorted.bam | head -1|cut -d " " -f1 > dedup
	samtools idxstats merged_bam/$i\_merged_dedup_sorted.bam | grep Ad5 |cut -f3 > ad5 
	paste name total uniq dedup ad5 > stat
	cat mapped_reads.txt stat > foo
	mv foo mapped_reads.txt
	rm total uniq dedup stat name ad5
done
sed -i '1isample\total_reads\tunique_reads\tdeduplicated\tAd5' mapped_reads.txt
