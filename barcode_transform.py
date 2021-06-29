import sys

filename = sys.argv[1]
# header = int(sys.argv[2])
#
# head = False
# if header:
#     head = True

count = 0
# first = '\t'.join(['index', 'fwd_bc', 'rev_bc', 'gene_num', 'UMI_num'])
# first = '\t'.join(['index', 'fwd_bc', 'rev_bc'])
# print first
with open(filename, 'r') as f:
    for l in f:
        # if head:
        #     head = False
        #     continue
        count += 1
        l = l.strip('\r\n')
        l = l.strip('\r')
        l = l.strip('\n')
        # l = l.split(',')
        # l = [i.strip('"') for i in l]
        rev = []
        for i in l:
            if i == 'A':
                rev.append('T')
            elif i == 'T':
                rev.append('A')
            elif i == 'G':
                rev.append('C')
            elif i == 'C':
                rev.append('G')
            else:
                print >> sys.stderr, "Something is wrong, try again", l
        out = ''.join(rev[::-1])
        # print '\t'.join([str(count), l[0], out, l[1], l[2]])
        print '\t'.join([str(count), l, out])
