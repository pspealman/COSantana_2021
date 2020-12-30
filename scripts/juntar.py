# -*- coding: utf-8 -*-
"""
juntar.py
Created on Wed Oct 21 22:06:48 2020

This is intended to replace join reads from Qiime which allowed read joining 
by pairs at less than 10 bp length. Additionally, if supplied, reads that can 
be resolved using vsearch will take precedence over read joining.

Standard usage:
    python juntar.py --fastq_1 /path_to_forward_read_file \
        --fastq_2 /path_to_reverse_read_file \
        --vsearch /path_to_vsearch_results_file \
        --output_file /path_for_output_files
        
v0.2:
    Added --minimum_overlap and --maximum_mismatch
        
"""
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-f1',"--fastq_1")
parser.add_argument('-f2',"--fastq_2")
parser.add_argument('-v',"--vsearch")
parser.add_argument('-min',"--minimum_overlap")
parser.add_argument('-err',"--maximum_mismatch")
parser.add_argument('-o',"--output_file")

args = parser.parse_args()

if args.minimum_overlap:
    minimum_overlap = int(args.minimum_overlap)
else:
    minimum_overlap = 4
    
if args.maximum_mismatch:
    maximum_mismatch = int(args.maximum_mismatch)
else:
    maximum_mismatch = 1

def reverse_compliment(oldstr):
    comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    newstr = ''
    isstr = oldstr[::-1]
    
    for each in isstr:
        newstr += (comp_dict[each])

    return(newstr)
    
def compare_two(f1_seq, f2_seq, f1_q, f2_q):
    global match_size
    
    for index in range(100, (minimum_overlap-1), -1):
        if f1_seq[(-1*index):] == f2_seq[:index]:
            
            f3_seq = f1_seq + f2_seq[index:]
            f3_q = f1_q + f2_q[index:]
            
            if len(f3_seq) != len(f3_q):
                print('seq v q error')
                print(len(f3_seq), len(f3_q))
                1/0
            
            if len(f3_seq) >= 250:
                if index not in match_size:
                    match_size[index] = 0
                
                match_size[index]+=1
            
                return(f3_seq, f3_q)               
            
    return(0, 0)
    
def slide_compare_two(f1_seq, f2_seq, f1_q, f2_q):
    hit=0
    miss=0
    
    for index in range(10):        
        if f1_seq[(-1*index)-1] == f2_seq[index]:
            hit += 1
        else:
            miss += 1
            
        if hit >= (minimum_overlap*2)-1:
            f3_seq = f1_seq + f2_seq[index:]
            f3_q = f1_q + f2_q[index:]
            
            if len(f3_seq) != len(f3_q):
                print('seq v q error')
                print(len(f3_seq), len(f3_q))
                1/0
            
            if len(f3_seq) >= 250:
                if index not in save_size:
                    save_size[index] = 0
                
                save_size[index]+=1
            
                return(f3_seq, f3_q) 
        
        if miss > maximum_mismatch:
            return(0,0)
            
    return(0, 0)

f1 = open(args.fastq_1)
f2 = open(args.fastq_2)
f3 = open(args.output_file, 'w')
f_log = open(args.output_file.split('.')[0] + '.log','w')
    
fi_dict = {}
fii_dict = {}

match_size = {}
total_hit = 0

save_size = {}
total_save = 0

v_set = set()
vmatch_set = set()

if args.vsearch:
    v1 = open(args.vsearch)
    
    ct = 0
    for line in v1:
        f3.write(line)
        ct+=1
        if ct == 1:
            read_id = line.split(' ')[0]
            v_set.add(read_id)
        if ct == 4:
            ct = 0 

    v1.close()

ct = 0
for line in f1:
    ct+=1
    if ct == 1:
        read_id = line.split(' ')[0]
        fi_dict[read_id]={'s':'', 'q':''}
    if ct == 2:
        fi_dict[read_id]['s'] = line.strip()
    if ct == 4:
        fi_dict[read_id]['q'] = line.strip()
        ct = 0  
f1.close()

ct = 0
for line in f2:
    ct+=1
    if ct == 1:
        read_id = line.split(' ')[0]
        fii_dict[read_id]={'s':'', 'q':''}
    if ct == 2:
        fii_dict[read_id]['s'] = reverse_compliment(line.strip())
    if ct == 4:
        line = line.strip()
        line = line[::-1]
        fii_dict[read_id]['q'] = line
        ct = 0  
f2.close()
    
for read_id in fi_dict:
    if (read_id in fii_dict):
        #if read_id =='@M01965:30:000000000-BMF6L:1:2114:15962:27429':
        #    1/0
        f1_seq = fi_dict[read_id]['s']
        f2_seq = fii_dict[read_id]['s']
        f1_q = fi_dict[read_id]['q']
        f2_q = fii_dict[read_id]['q']        
        
        f3_seq, f3_q = compare_two(f1_seq, f2_seq, f1_q, f2_q)
        
        if f3_seq != 0:
            if (read_id not in v_set):
                total_hit+=1
                #read_id = read_id + ' 1:N:0:GGACTCCT+GCGTAAGA'
                outline = ('{read_id}\n{f3_seq}\n+\n{f3_q}\n').format(read_id=read_id, f3_seq=f3_seq, f3_q =f3_q)
                f3.write(outline)
                
            if (read_id in v_set):
                vmatch_set.add(read_id)
            
        else:
            f3_seq, f3_q = slide_compare_two(f1_seq, f2_seq, f1_q, f2_q)
            
            if f3_seq != 0:
                total_save+=1
                #read_id = read_id + ' 1:N:0:GGACTCCT+GCGTAAGA'
                outline = ('{read_id}\n{f3_seq}\n+\n{f3_q}\n').format(read_id=read_id, f3_seq=f3_seq, f3_q =f3_q)
                f3.write(outline)
                
                if (read_id in v_set):
                    vmatch_set.add(read_id)

f3.close()

outline = ('\tTotal Reads:\t{total}\n\tTotal import from Vsearch:\t{vsearch}\t({frac}%)\n').format(total = len(fi_dict), vsearch = len(v_set), frac = len(v_set)/len(fi_dict))
print(outline)
f_log.write(outline)

outline = ('\tTotal Perfect Join Match:\t{join}\t({frac}%)\n').format(join = total_hit, frac = total_hit/len(fi_dict))
print(outline)
f_log.write(outline)

outline = ('\tTotal Saved Join:\t{save}\t({frac}%)\n').format(save = total_save, frac = total_save/len(fi_dict))
print(outline)
f_log.write(outline)

f_log.close()
