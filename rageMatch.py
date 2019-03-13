import sys, os, re, gzip
import argparse

'''
    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2018

    Direct matching and trimming of fastq files

    sge:

    source ~/venv2714/bin/activate
    time python rageMatch.py -q ramos.fastq.gz -b Ramos_master_barcodes.tsv -n ${SGE_TASK_ID} -c 4000 -t 10 -o ${TMPDIR}/ramos -l 0 -f 250 -s
    cp ${TMPDIR}/ramos/* ./output

    Command:

    LINES=$(wc -l Ramos_master_barcodes.tsv | awk '{print $1}')

    CMD="qsub -cwd -V -pe smp 1 -N R_DM -S /bin/bash -t 1-$LINES -l mem_requested=8G,h_vmem=8G,tmp_requested=50G ../rageMatch.sge"
    echo $CMD && $CMD

'''


def bc_pull(args):
    bc_num = args.barcodeNumber - 1
    bc_list = []
    with open(args.barcodes , 'r') as bc:
        for line in bc:
            line = line.strip('\n')
            l_list = line.split('\t')
            bc_list.append(l_list)
    return bc_list[bc_num]

def main():

    parser = argparse.ArgumentParser(description="rageMatch - RageSeq direct matching and trimming")
    parser.add_argument("-q", "--fastq",
                        help="fastq file - gzipped")
    parser.add_argument("-o", "--output",
                        help="barcode fastq output prefix")
    parser.add_argument("-b", "--barcodes",
                        help="barcode file - built from transform_bc.py")
    parser.add_argument("-n", "--barcodeNumber", type=int,
                        help="The barcode line number from the barcode file 1 indexed - used for HPC array jobs or a single barcode")
    parser.add_argument("-c", "--chunks", type=int,
                        help="Number of reads to store in memory before processing")
    parser.add_argument("-l", "--length", type=int,
                        help="Number of nt to search at start and end")
    parser.add_argument("-f", "--filter", type=int,
                        help="Min length of read required to be used")
    parser.add_argument("-t", "--trim", type=int,
                        help="Trim barcodes and adapters from reads - number of bases")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Engage higher output verbosity")
    parser.add_argument("-s", "--stats", action="store_true",
                        help="dump match stats into header")
    args = parser.parse_args()

    out = args.output +'.' + str(args.barcodeNumber) + '.fastq'
    bcs = bc_pull(args)
    fq_list = []
    fq_tmp = []
    chunk = 0
    with gzip.open(args.fastq, 'rb') as fq_file:
        c = 0
        for line in fq_file:
            c += 1
            fq_tmp.append(line.strip())
            if c == 4:
                chunk += 1
                fq_list.append(fq_tmp)
                fq_tmp = []
                c = 0
                if chunk >= args.chunks:
                    do_batch(args, fq_list, bcs, out)
                    fq_list = []
                    chunk = 0

    do_batch(args, fq_list, bcs, out)

def do_batch(args, faq, bc, out):
    with open(out, 'a') as o1:
        for fq in faq:
            if len(fq[1]) < args.filter:
                continue
            if args.length == 0:
                if len(fq[1]) > 400:
                    length = int(len(fq[1]) / 2.0)
                else:
                    length = 200
            k = re.search(bc[2], fq[1][:length])
            n = re.search(bc[3], fq[1][len(fq[1]) - length:])
            if k:
                # +13 to give a buffer of 3 for the UMI into the poly TAil
                if args.trim:
                    trim = k.end() + (args.trim + 3)
                    if args.stats:
                        fq[0] = '{} fwd_bc={} fwd_UMI={} cellNum={} start={} end={}'.format(fq[0], bc[2], fq[1][k.end():k.end()+args.trim], str(args.barcodeNumber), str(k.start()), str(k.end()))
                    else:
                        fq[0] = '{} fwd_bc={} fwd_UMI={} cellNum={}'.format(fq[0], bc[2], fq[1][k.end():k.end()+args.trim], str(args.barcodeNumber))
                    fq[1] = fq[1][trim:]
                    fq[3] = fq[3][trim:]
                else:
                    fq[0] = '{} fwd_bc={} cellNum={}'.format(fq[0], bc[2], str(args.barcodeNumber))
                o1.write('\n'.join(fq))
                o1.write('\n')
            elif n:
                # only -12 because of where it stops so goes one further, still 3 buffer
                if args.trim:
                    trim = n.start() - (args.trim + 2)
                    trim = len(fq[1]) - (length - trim)
                    if args.stats:
                        fq[0] = '{} fwd_bc={} fwd_UMI={} cellNum={} start={} end={}'.format(fq[0], bc[3], fq[1][trim+2:trim+args.trim+2], str(args.barcodeNumber), str(len(fq[1]) - length + n.start()), str(len(fq[1]) - length + n.end()))
                    else:
                        fq[0] = '{} rev_bc={} rev_UMI={} cellNum={}'.format(fq[0], bc[3], fq[1][trim+2:trim+args.trim+2], str(args.barcodeNumber))
                    fq[1] = fq[1][:trim]
                    fq[3] = fq[3][:trim]
                else:
                    fq[0] = '{} rev_bc={} cellNum={}'.format(fq[0], bc[3], str(args.barcodeNumber))
                o1.write('\n'.join(fq))
                o1.write('\n')


if __name__ == '__main__':
    main()
