#!/bin/bash
#Copyright 2019 Martin A. Smith, m.smith[at]garvan.org.au

READS=$( ls lymph/cells/*fastq | sort -h -t "." -k 2 | sed -n ${SGE_TASK_ID}p  )
T_READS=${READS##*/}
SUMMARY=lymph/seq_sum.txt.gz
INDEX=lymph/lymph.index.gz
THREADS=4

if [[ ! -d logs ]] ; then mkdir logs ;  fi
#if [[ ! -d bams ]] ; then mkdir bams ;  fi
if [[ ! -d fasta ]] ; then mkdir fasta ; fi

if [[ ! -e ${SUMMARY} ]]
then
  >&2 echo -e "[ERROR] Sequencing Summary file not found: "${SUMMARY}"\nAborting"
  exit 1
fi

if [[ ! -e ${INDEX} ]]
then
  >&2 echo -e "[ERROR] Fast5 index file not found: "${INDEX}"\nAborting"
  exit 1
fi


>&2 echo "[*] Launching contig assembly for "$READS
canu  -d ${TMPDIR}/canu \
  -p test \
  -useGrid=false \
  -Overlapper=minimap \
  -batMemory=28 \
  -minReadLength=200 \
  -minOverlapLength=100 \
  -genomeSize=15k \
  -maxThreads=${THREADS} \
  -minThreads=${THREADS} \
  -nanopore-raw ${READS} 2> logs/canu.${SGE_TASK_ID}.log
# -stopOnReadQuality=false \

DRAFT=${TMPDIR}/canu/test.contigs.fasta
if [[ ! -s ${DRAFT} ]] ; then
  echo "[*] No contigs assembled... exiting"
  qstat -j ${JOB_ID} | grep usage
  exit 0
fi

cp $DRAFT fasta/${T_READS%*.fastq}_contigs_canu.${SGE_TASK_ID}.fa
echo "[*] Launching consensus calculation"
for ROUND in {01..04}; do
  echo "Running round ${ROUND} consensus..."
  READS2TIGS=${TMPDIR}/reads2contigs_${ROUND}.paf
  NEWDRAFT=${TMPDIR}/racon_${ROUND}.fasta
  minimap2 -k 15 -t ${THREADS} ${DRAFT} ${READS} > ${READS2TIGS}
  racon -w 200 -m 8 -x -1 -g -4  -t ${THREADS} -q -1 ${READS} ${READS2TIGS} ${DRAFT} > ${NEWDRAFT}
  #rm ${DRAFT}
  DRAFT=${NEWDRAFT}
done 2> logs/consensus.${SGE_TASK_ID}.log
cp $DRAFT fasta/${T_READS%*.fastq}_contigs_racon.${SGE_TASK_ID}.fa

# echo "[*] mapping reads to consensus assembly"
# BAM=$(pwd)/bams/${T_READS%*.fastq}_map2contig.${SGE_TASK_ID}.bam
# ( minimap2 -k 15 -a -t ${THREADS}  ${DRAFT} ${READS} | samtools view -b - | samtools sort -@ 3 -o ${BAM} ) 2>&1 >> logs/minimap.${SGE_TASK_ID}.log
# samtools index ${BAM} 2>&1 > logs/samtools.${SGE_TASK_ID}.log
mv ${DRAFT} fasta/${T_READS%*.fastq}_contigs_racon.${SGE_TASK_ID}.fa
rm ${TMPDIR}/racon* ${TMPDIR}/reads2contigs*

echo "All done! Usage: "
qstat -j ${JOB_ID} | grep usage
