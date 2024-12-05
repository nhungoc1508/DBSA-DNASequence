import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Compare query results with and without SP-GiST index")
parser.add_argument('--dir', '-d', help='Directory', required=True, type=str)
args = parser.parse_args()

def compare(s, t):
    t = list(t)   # make a mutable copy
    try:
        for elem in s:
            t.remove(elem)
    except ValueError:
        print(f'In S not in T: {elem}')
        return False
    return not t

def main():
    files = []
    for f in os.listdir(f'./{args.dir}'):
        files.append(f)
    if len(files) != 2:
        print('There must be exactly two files for comparison')
        sys.exit(1)
    seqscan = []
    idxscan = []
    for file in files:
        with open(f'./{args.dir}/{file}', 'r') as in_file:
            lines = in_file.readlines()
            for i, line in enumerate(lines):
                if 'kmer' in line and '------' in lines[i+1]:
                    idx = i + 2
                    break
            for i in range(idx, len(lines)):
                if 'row' not in lines[i]:
                    line = lines[i].strip()
                    if len(line) > 0:
                        if 'seq' in file:
                            seqscan.append(line)
                        else:
                            idxscan.append(line)
    
    if (compare(seqscan, idxscan)):
        print(f'Compare seq - idx: Results are identical, number of rows returned: {len(seqscan)}')
    else:
        print('Results are not the same')

    if (compare(idxscan, seqscan)):
        print(f'Compare idx - seq: Results are identical, number of rows returned: {len(idxscan)}')
    else:
        print('Results are not the same')

if __name__ == '__main__':
    main()