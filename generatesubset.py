import sys
import re
import numpy as np

gtf = './inputs/gencode.v32.primary_assembly.annotation.chr22.gtf'
fasta = './inputs/GRCh38.primary_assembly.genome.chr22.fa' #sys.argv[2]
bam = './Aligned.out.sam' #sys.arv[3]

class BamRanges:
    def getsimpleflanks(self, bam):
        bamffile = open(bam,'r')
        maps = {}
        for line in bamffile:
            if line[0] == '@':
                continue
            tokens = line.split('\t')
            chrom = tokens[2]
            pos = int(tokens[3])
            length = sum(map(int, re.findall(r'\d+', tokens[5])))

            if chrom not in maps:
                maps[chrom] = []

            maps[chrom].append( (pos, pos + length) )
        return(maps)

    def flatten(self, flank_map):
        maps = {}
        for chrom in flank_map:
            ranges = {}
            for beg,end in flank_map[chrom]:
                ranges[beg] = True
                ranges[end] = True

            data = sorted(ranges.keys())
            maps[chrom] = data
        return(maps)

    def __init__(self, bam, thresh):
        self.ranges = self.flatten(self.getsimpleflanks(bam))
        self.clustered = self.clusterbreaks_threshold(self.ranges, thresh)

    # threshold method
    def clusterbreaks_threshold(self, arrs, thresh = 100):
        isles = {}
        for chrom in arrs:
            array = arrs[chrom]
            islands = [[array[0]]]
            gaps = []
            i = 1
            while i < len(array):
                num = array[i]
                last = islands[-1][-1]
                diff = num - last

                if diff <= thresh:
                    islands[-1].append(num)
                else:
                    islands[-1].append(last + 200)
                    gaps.append(diff)
                    islands.append([num])
                i += 1
            print(len(islands), ':', ' '.join(map(str,sorted([len(x) for x in islands]))))
            print('\n'.join(map(str, sorted(gaps))))
            isles[chrom] = islands
        return(isles)


    # GMM method -- fix there
    from sklearn.mixture import GaussianMixture
    def clusterbreaks_GM(array):
        X = np.array(array).reshape(-1,1)
        model = GaussianMixture(n_components=40, covariance_type='spherical')
        model = model.fit(X)
        # compute the AIC and the BIC
        #AIC = [m.aic(X) for m in models]
        #BIC = [m.bic(X) for m in models]



def fastasubset(fasta, chrom_ranges):
    fastafile = open(fasta, 'r')
    fastaout = open(fasta + '.subsetted.fa', 'w')
    active_ranges = False
    current_index = 0
    thresh = 200
    for line in fastafile:
        line = line.splitlines()[0]
        if line[0] == '>':
            # Header
            current_index = 0
            chrom = line[1:].split(' ')[0]
            if chrom in chrom_ranges:
                active_ranges = chrom_ranges[chrom]
                print(line, "using")
                print(line, file=fastaout)
            else:
                active_ranges = False
                print(line, "skip")
        else:
            # Region
            if not active_ranges:
                continue

            # active region
            for cluster in active_ranges:
                # use extremes of the clusters
                beg, end = (cluster[0], cluster[-1])
                if (beg-thresh) <= current_index <= (end+thresh):
                    print(line, file=fastaout)

            current_index += len(line)


def fastasubset_exact(fasta, chrom_ranges):
    '''Method that extracts exact positions instead of the
    per-line positions given previously'''
    fastafile = open(fasta, 'r')
    fastaout = open(fasta + '.subsetted.exact.fa', 'w')

    # Reset these per crhom
    active_ranges = False
    active_range_data = '' # concat all lines of same chrom

    def processactivedata():
        if active_range_data != '':
            extracted = []
            for rang in active_ranges:
                beg = rang[0]
                end = rang[-1]

                sl = slice(beg + 1,end + 1,1) # 1 -based numbering
                #print(beg,end,1)
                sub = active_range_data[sl]
                extracted.append(sub)

            outstring = ''.join(extracted)
            print("out sequence", len(outstring), "characters")
            print(outstring, file=fastaout)


    for line in fastafile:
        line = line.splitlines()[0]
        if line[0] == '>':
            # Header
            current_index = 0
            chrom = line[1:].split(' ')[0]
            if chrom in chrom_ranges:
                active_ranges = chrom_ranges[chrom]
                print(line, "using")
                print(line, file=fastaout)
            else:
                # If activate_range_data is full or we are EOF then we have
                # string from previous that needs to be processed
                processactivedata()

                # reset
                active_ranges = False
                active_range_data = ''
                print(line, "skip")
        else:
            # Region
            if not active_ranges:
                continue

            # active region -- concat all into supermassive string
            active_range_data += line

    processactivedata()



def getpositionmap(chrom_ranges):
    '''Converts the chromosome ranges to their
    new positions given by the subsetted fasta which
    are now shifted'''
    posmap = {}
    for chrom in chrom_ranges:
        # In a chrom
        intervals = [0]
        actual = [chrom_ranges[chrom][0][0]]
        for rang in chrom_ranges[chrom]:
            beg = rang[0]
            end = rang[-1]
            interval = end - beg
            intervals.append(intervals[-1] + interval)
            actual.append(end)
        posmap[chrom] = {'interval' : intervals, 'actual' : actual}
    return(posmap)



def gtfsubset_exact(gtf, posmap):
    with open(gtf, 'r') as gtffile:
        gtfout = open(gtf + '.subsetted.exact.gtf', 'w')
        
        for line in gtffile:
            line = line.splitlines()[0]
            if line[0]=='#':
                print(line, file=gtfout)
                continue
            
            tokens = line.split('\t')
            chrom, beg, end = tokens[0], int(tokens[3]), int(tokens[4])
            if chrom in posmap:
                fasta_positions = posmap[chrom]['interval']
                actual_positions = posmap[chrom]['actual']

                # Find the interval that beg and end are in
                i = 0
                while i + 1 < len(actual_positions):
                    aposi = actual_positions[i]
                    anext = actual_positions[i+1]
                    
                    found = 0
                    if aposi <= beg <= anext:
                        found += 1
                    if aposi <= end <= anext:
                        found += 1

                    # the whole feature is in the interval
                    # let's just take only these features, otherwise it gets complex...
                    if found == 2:
                        # Get updated interval
                        offset = fasta_positions[i]
                        new_beg = (beg - aposi) + offset
                        new_end = (end - aposi) + offset
                        
                        tokens[3] = str(new_beg)
                        tokens[4] = str(new_end)
                        
                        print('\t'.join(tokens), file=gtfout)

                    i += 1



# {'chr22': (10953092, 50798952)}
# We choose 150000 to keep the file size less than 1M
chrom_ranges = BamRanges(bam, 150000).clustered
# Generate the new fasta
fastasubset_exact(fasta, chrom_ranges)
# Generate the new GTF
