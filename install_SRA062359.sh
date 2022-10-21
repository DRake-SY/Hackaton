for i in "636533" "636531" "636532" "628589" "628588" "628587" "628586" "628585" "628584" "628583" "628582";
do
prefetch SRR${i};
fastq-dump --gzip --split-files SRR${i}/SRR${i}.sra;
done;
