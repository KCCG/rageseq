#!/bin/bash
#Copyright 2019 Martin A. Smith, m.smith[at]garvan.org.au

READS=$( ls lymph/cells/*fastq | sort -h -t "." -k 2 | sed -n ${SGE_TASK_ID}p  )
T_READS=${READS##*/}
SUMMARY=lymph/seq_sum.txt.gz
INDEX=lymph/lymph.index.gz
THREADS=1

if [[ ! -d logs ]] ; then mkdir logs ;  fi
if [[ ! -d bams ]] ; then mkdir bams ;  fi
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

if [[ ! -e fasta/${T_READS%*.fastq}_contigs_racon.${SGE_TASK_ID}.fa  ]]
then
  >&2 echo "[*] No consensus for fasta/${T_READS%*.fastq}_contigs_racon.${SGE_TASK_ID}.fa"
  exit 1
fi

if [[ -s fasta/${T_READS%*.fastq}_contigs_nanopolish.${SGE_TASK_ID}.fa  ]]
then
  >&2 echo "[*] Already polished it off. ABORT!"
  exit 1
fi

CONTIG=fasta/${T_READS%*.fastq}_contigs_racon.${SGE_TASK_ID}.fa

echo "[*] mapping reads to consensus assembly"
BAM=$(pwd)/bams/${T_READS%*.fastq}_map2contig.${SGE_TASK_ID}.bam
( minimap2 -k 15 -a -t ${THREADS}  ${CONTIG} ${READS} | samtools view -b - | samtools sort -@ 3 -o ${BAM} ) 2>&1 >> logs/minimap.${SGE_TASK_ID}.log
samtools index ${BAM} 2>&1 > logs/samtools.${SGE_TASK_ID}.log

# Check if fast5.tar exists
if [[ ! -e ${T_READS%*.fastq}_fast5s.tar  ]]
then
  echo "[*] Extracting fast5s"
  if [[ ! -e ${TMPDIR}/fast5 ]]; then mkdir ${TMPDIR}/fast5 ; fi
  CMD="python fast5_fetcher.py -q ${READS} -s ${SUMMARY} -i ${INDEX} -o ${TMPDIR}/fast5"
  echo $CMD && $CMD
#  cd ${TMPDIR}
#  tar -cf ${T_READS%*.fastq}_fast5s.tar fast5/
#  mv ${T_READS%*.fastq}_fast5s.tar ${SGE_O_WORKDIR}/
#  cd ${SGE_O_WORKDIR}
#else
#  cd ${TMPDIR}
#  tar -xf ${SGE_O_WORKDIR}/${T_READS%*.fastq}_fast5s.tar
#  cd ${SGE_O_WORKDIR}
fi

echo "[*] Launching nanopolish"
CMD="nanopolish index -d ${TMPDIR}/fast5 ${READS}"
echo $CMD && $CMD 2>&1 > logs/nanopolish.${SGE_TASK_ID}.log

# $cat ramos_cell_1001_contigs.2.fa.fai
# tig00000001    1310    46    1310    1311
# tig00000003    968    1401    968    969
# tig00000005    415    2414    415    416

CMD="samtools faidx ${CONTIG}"
echo $CMD && $CMD

if [[ ! -e ${CONTIG}.fai ]]
then
  echo "[ERROR] Fasta indexing failed. Aborting"
  exit 1
fi

cat  ${CONTIG}.fai | while read LINE
do
  CONT_ID=$( awk '{print $1}' <<< $LINE )
  echo "[**] Polishing contig "${CONT_ID}
  CONT_LEN=$( awk '{print $2}' <<< $LINE )
  CONT_LEN=$(( $CONT_LEN + 1 ))
  CMD="nanopolish variants -t ${THREADS}"
  CMD=$CMD" -w "$CONT_ID":1-"${CONT_LEN}
  CMD=$CMD" --consensus=${TMPDIR}/${CONT_ID}-nanop.fa"
  CMD=$CMD" --reads ${READS}"
  CMD=$CMD" --bam ${BAM}"
  CMD=$CMD" --genome ${CONTIG}"
  echo $CMD && $CMD 2>&1 >> logs/nanopolish.${SGE_TASK_ID}.log
done
cat ${TMPDIR}/*-nanop.*fa > ${CONTIG%_*}_nanopolish.${SGE_TASK_ID}.fa

rm ${BAM} ${BAM}.bai ${READS}.index* ${CONTIG}.fai

echo "All done! Usage: "

qstat -j ${JOB_ID} | grep usage





mv ${DRAFT} fasta/${T_READS%*.fastq}_contigs_racon.${SGE_TASK_ID}.fa
rm ${TMPDIR}/racon* ${TMPDIR}/reads2contigs*

echo "All done! Usage: "
qstat -j ${JOB_ID} | grep usage
