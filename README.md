# rageseq
RAGE-seq scripts
#### Barcode extraction


#### Demultiplexing

##### rageMatch.py

Used to read the barcode whitelist, search the total fastq, bin reads on barcode matches, and trim to poly TAil.

Used on an HPC environment with Sun Grid Engine


sge:

    source ~/venv2714/bin/activate
    time python rageMatch.py -q ramos.fastq.gz -b Ramos_master_barcodes.tsv -n ${SGE_TASK_ID} -c 4000 -t 10 -o ${TMPDIR}/ramos -l 0 -f 250 -s
    cp ${TMPDIR}/ramos/* ./ramos/output

Command:

    LINES=$(wc -l Ramos_master_barcodes.tsv | awk '{print $1}')

    CMD="qsub -cwd -V -pe smp 1 -N R_DM -S /bin/bash -t 1-$LINES -l mem_requested=8G,h_vmem=8G,tmp_requested=50G ../rageMatch.sge"
    echo $CMD && $CMD



```
usage: rageMatch.py [-h] [-q FASTQ] [-o OUTPUT] [-b BARCODES]
                    [-n BARCODENUMBER] [-c CHUNKS] [-l LENGTH] [-f FILTER]
                    [-t TRIM] [-v] [-s]

rageMatch - RageSeq direct matching and trimming

optional arguments:
  -h, --help            show this help message and exit
  -q FASTQ, --fastq FASTQ
                        fastq file - gzipped
  -o OUTPUT, --output OUTPUT
                        barcode fastq output prefix
  -b BARCODES, --barcodes BARCODES
                        barcode file - built from transform_bc.py
  -n BARCODENUMBER, --barcodeNumber BARCODENUMBER
                        The barcode line number from the barcode file 1
                        indexed - used for HPC array jobs or a single barcode
  -c CHUNKS, --chunks CHUNKS
                        Number of reads to store in memory before processing
  -l LENGTH, --length LENGTH
                        Number of nt to search at start and end
  -f FILTER, --filter FILTER
                        Min length of read required to be used
  -t TRIM, --trim TRIM  Trim barcodes and adapters from reads - number of
                        bases
  -v, --verbose         Engage higher output verbosity
  -s, --stats           dump match stats into header
```    



#### Assebmly



See our preprint:

https://www.biorxiv.org/content/early/2018/09/24/424945
