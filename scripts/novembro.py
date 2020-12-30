# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 09:29:45 2019

Novembro - for identifying taxa enrichments from Qiime generated dat

python novembro.py -f feature-table.biom.txt -t taxonomy.tsv -s silva -o taxa_counts.tab

NB: feature-table.biom.txt should be generated from converting the qiime feature_table.biom file to a tsv
#   biom convert -i feature-table.biom -o feature-table.biom.txt --to-tsv

Version: Public 1.0 (I fill my day with hope and face it with joy.)
Version: Public 1.1 (Piano Harmony)
    _x_ added Kruskal-Wallis, --kruskal_wallis
    _x_ kludge taxa_raw_dict iteration to handle non-triplicate sets
Version: Public 1.2 (Diamond Retiree)
    _x_ added unique_function
Version: Public 1.3 (Formula Choice)
    _x_ added restructuring for 4 samples with 3 replicates   
Version: Public 1.4 (Art Charity)
    _x_ make coherent for BGS
    _x_ added arguments for 'pct_effect_size, pval_threshold'
    _x_ added binomial exact test for abundances
    
@author: ps163@nyu.edu
"""

import numpy as np
import scipy.stats as stats

import plotly.graph_objects as go
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f',"--feature_table")
parser.add_argument('-t',"--taxonomy")
parser.add_argument('-s',"--taxa_source")
parser.add_argument('-o',"--output_file")

parser.add_argument('-u',"--unique_sets")

parser.add_argument('-pct',"--pct_effect_size")
parser.add_argument('-pval',"--pval_threshold")
parser.add_argument('-kw',"--kruskal_wallis", action='store_true')

args = parser.parse_args()

if args.pval_threshold:
    pval_threshold = float(args.pval_threshold)
else:
    pval_threshold = 0.05
    
if args.pct_effect_size:
    pct_effect_size = float(args.pct_effect_size)
else:
    pct_effect_size = 0.05
    
if args.kruskal_wallis:
    stats_runmode = 'kruskal_wallis'
else:
    stats_runmode = 'chi2'
    
feature_table_name = args.feature_table 
taxa_file_name = args.taxonomy

if args.taxa_source:
    taxa_source = args.taxa_source
else:
    taxa_source = 'silva'

if taxa_source.lower() == 'silva' or taxa_source.lower() == 's':
    prefixe = ['D_0__','D_1__','D_2__','D_3__','D_4__','D_5__','D_6__','D_7__','D_8__','D_9__','D_10__','D_11__','D_12__','D_13__','D_14__']
   
if 'green' in taxa_source.lower() or taxa_source.lower() == 'gg':
    prefixe = ['k__','p__','c__','o__','f__','g__','s__']
    
rank_order = ['species','genus','family','order','class','phylum','kingdom']
convert_taxa_to_rank = {'kingdom':0, 'phylum':1, 'class':2, 'order':3, 'family':4, 'genus':5, 'species': 6}

#filter_previous = []
pass_dict = {'criteria':0, 'pval':0, 'p_pval':0, 'pass_set':0, 'figure_dict':0}

simplified_enrichment = {}

def increment_pass_dict(p_1v2, p_1v3=1, p_2v3=1):
    xlist = [p_1v2, p_1v3, p_2v3]
    
    for each_p in xlist:
        if each_p <= 0.05:
            return(1)
    return(0)
    
#def evaluate_the_mwu(p_val, w_val, testname, max_effect_size, min_effect_size, pval_threshold, pass_set):
#    if (p_val <= pval_threshold) and (w_val >= max_effect_size or w_val <= min_effect_size):
#        pass_set.append(testname)
#    return(pass_set)
    
def evaluate_the_criteria(p_val, w_val, testname, max_effect_size, min_effect_size, pval_threshold, pass_set):
    if (p_val <= pval_threshold) and (w_val >= max_effect_size or w_val <= min_effect_size):
        pass_set.append(testname)
    return(pass_set)
    
def obs_counter(list_obs_ct):
    obs = 0
    
    if max(list_obs_ct) > 100:
        return(True)
        
    for obs_ct in list_obs_ct:
        if obs_ct > 100:
            obs+=1
            
    if obs > 1:
        return(True)
    else:
        return(False)

def find_correction_value(otu_counts):   
    sub1, sub2, sub3 = 0, 0, 0
    int1, int2, int3 = 0, 0, 0
    sup1, sup2, sup3 = 0, 0, 0
    
    for otu in otu_counts:
       osub1, osub2, osub3, oint1, oint2, oint3, osup1, osup2, osup3 = otu_counts[otu]
       
       sub1 += osub1
       sub2 += osub2
       sub3 += osub3
       
       int1 += oint1
       int2 += oint2
       int3 += oint3
       
       sup1 += osup1
       sup2 += osup2
       sup3 += osup3
                   
    global_min = min(min([sub1, sub2, sub3]),
                     min([int1, int2, int3]),
                     min([sup1, sup2, sup3])
                     )
    
    #Define number of observations per site for the purposes of downsampling
    ct_cor_1, ct_cor_2, ct_cor_3 = global_min/sub1, global_min/sub2, global_min/sub3
    g_cor_1, g_cor_2, g_cor_3 = global_min/int1, global_min/int2, global_min/int3  
    ol_cor_1, ol_cor_2, ol_cor_3 = global_min/sup1, global_min/sup2, global_min/sup3
    
    return(ct_cor_1, ct_cor_2, ct_cor_3, g_cor_1, g_cor_2, g_cor_3, ol_cor_1, ol_cor_2, ol_cor_3)

def build_otu_counts(feature_table_name):
#this is to build the 
    otu_file = open(feature_table_name)
    
    otu_counts = {}
    
    for line in otu_file:
        if line[0]!='#' and 'Feature ID' not in line:
            #OTU ID	sub1	sub2	sub3	int1	int2	int3	sup1	sup2	sup3	S1	S2	S3
            line = line.strip()
            otu = line.split('\t')[0].strip()# 
            #
            sub1 = float(line.split('\t')[1])
            sub2 = float(line.split('\t')[2])  
            sub3 = float(line.split('\t')[3])
            #
            int1 = float(line.split('\t')[4])
            int2 = float(line.split('\t')[5])  
            int3 = float(line.split('\t')[6])
            #
            sup1 = float(line.split('\t')[7])
            sup2 = float(line.split('\t')[8])  
            sup3 = float(line.split('\t')[9])            
            otu_counts[otu] = [sub1, sub2, sub3, int1, int2, int3, sup1, sup2, sup3]
    #        
    otu_file.close()
    return(otu_counts)
    
def criteria(p_set, s_set, v_set, p_1v2, p_1v3, p_2v3, pval, taxa, pct_effect_size=0.05, pval_threshold=0.05):
    pass_set = []
    
    global pass_dict
        
    pass_dict['criteria']+=1
    #print(taxa)
    if (sum(p_set)+ sum(s_set) + sum(v_set)) >= 100:
        if (pval <= pval_threshold):
            pass_dict['pval']+=1
            max_effect_size = (1+pct_effect_size)
            min_effect_size = (1-pct_effect_size)
            
            p_mean = np.mean(p_set)
            s_mean = np.mean(s_set)
            v_mean = np.mean(v_set)
                    
            if p_mean == 0:
                p_mean = 1
            if s_mean == 0:
                s_mean = 1
            if v_mean == 0:
                v_mean = 1  
                              
            w_1v2 = (p_mean/s_mean)
            w_1v3 = (p_mean/v_mean)
            w_2v3 = (s_mean/v_mean)
            
            #log_fold_diff = [w_1v2, w_1v3, w_1v4, w_2v3, w_2v4, w_3v4]
            
            pass_dict['p_pval'] += increment_pass_dict(p_1v2, p_1v3, p_2v3)
            #print((p_1v2, p_1v3, p_2v3))
            #print(w_1v2, w_1v3, w_2v3)
            #
            pass_set = evaluate_the_criteria(p_1v2, w_1v2, '1v2', max_effect_size, min_effect_size, pval_threshold, pass_set)
            print(pass_set)
            pass_set = evaluate_the_criteria(p_1v3, w_1v3, '1v3', max_effect_size, min_effect_size, pval_threshold, pass_set)
            print(pass_set)
            pass_set = evaluate_the_criteria(p_2v3, w_2v3, '2v3', max_effect_size, min_effect_size, pval_threshold, pass_set)
            print(pass_set)
                    
        if len(pass_set) >= 2:
            return(True)            

    return(False)
        
def test_max(set_list, max_obs):
    for each_set in set_list:
        if max(each_set) >= max_obs:
            max_obs = max(each_set)
    return(max_obs)

def return_log10(each_set, fraction_correction=False):
    new_set = []
    
    for each_obs in each_set:
        if fraction_correction:
            if each_obs < 1:
                each_obs = 0
            else:
                each_obs = np.log10(each_obs)
        else:
            if each_obs == 0:
                each_obs = 0
            else:
                each_obs = np.log10(each_obs)

        new_set.append(each_obs)
                
    return(new_set)
    
def return_deets(x_array, y_array):
    return(x_array, np.mean(x_array), np.std(x_array), y_array, obs_counter(y_array))
    
def run_mwu(x, y, z):
    pval_x_y, pval_x_z, pval_y_z = 1, 1, 1
    
    if (x != y):
        if ((sum(x) > 20) and sum(y) > 20):
            _ps, pval_x_y = stats.mannwhitneyu(x, y)
    if (x != z):
        if ((sum(x) > 20) and sum(z) > 20):
            _pv, pval_x_z = stats.mannwhitneyu(x, z)
    if (x != z):
        if ((sum(y) > 20) and sum(z) > 20):
            _sv, pval_y_z = stats.mannwhitneyu(y, z)
    
    return(pval_x_y, pval_x_z, pval_y_z)
    
def run_kruskal(x_set, y_set):
    if sum(x_set) > 30 or sum(y_set) > 30:
        _w, p_xvy = stats.kruskal(x_set, y_set)
        return(p_xvy)
    else:
        return(1)
    
def run_kruskal_3way(x_set, y_set, z_set):
    if sum(x_set) > 30 or sum(y_set) > 30 or sum(z_set) > 30:
        _w, p_xvy = stats.kruskal(x_set, y_set, z_set)
        return(p_xvy)
    else:
        return(1)
        
def run_chi2(x_num, x_den, y_num, y_den):
    if (x_num) > 5 or (y_num) > 5:
        obs = np.array([[max(x_num,1), x_den], [max(y_num,1), y_den]])
        chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
        return(pval)
    else:
        return(1)
        
def mod_null_set(isset):
    if len(isset) <1:
        isset.append(0)
    return(isset)
    
def simplify_enrichment(taxa, p12, p13, p23, p_med, s_med, v_med):
    global simplified_enrichment
    
    print(taxa, p12, p13, p23, p_med, s_med, v_med)
    simplified_enrichment[taxa] = 'complex'
    
    if max(p12,p13,p23)<=0.05:
        if p_med > s_med and p_med > v_med:
            simplified_enrichment[taxa]='Sub_high'
        if s_med > p_med and s_med > v_med:
            simplified_enrichment[taxa]='Inter_high'
        if v_med > p_med and v_med > s_med:
            simplified_enrichment[taxa]='Supra_high'        
    else:
        if p12 <= 0.05 and p13 <= 0.05:
            if p_med > s_med and p_med > v_med:
                simplified_enrichment[taxa]='Sub_high'
            if p_med < s_med and p_med < v_med:
                simplified_enrichment[taxa]='Sub_low'
    
        if p12 <= 0.05 and p23 <= 0.05:
            if s_med > p_med and s_med > v_med:
                simplified_enrichment[taxa]='Inter_high'
            if s_med < p_med and s_med < v_med:
                simplified_enrichment[taxa]='Inter_low'           

        if p13 <= 0.05 and p23 <= 0.05:
            if v_med > p_med and v_med > s_med:
                simplified_enrichment[taxa]='Supra_high'
            if v_med < p_med and v_med < s_med:
                simplified_enrichment[taxa]='Supra_low'
        
    return(simplified_enrichment[taxa])
    
    
def plot_top10_taxa(all_taxa_dict, prefix_name, pct_threshold=0.01):
    #taxa_dict[taxa][0] += sum(pristine)
    #taxa_dict['total'][0] += sum(pristine)
    
    logfile_name = ('taxonomic_abundance_{}_{}_{}.log').format(prefix_name, pct_threshold, taxa_cutoff_name)
    logfile = open(logfile_name, 'w')
    header = ('taxa\tsub_log10_abundance\tint_log10_abundance\tsup_log10_abundance\n')
    logfile.write(header)
    
    specific_taxa_dict = {'taxa':[],'CT':[], 'G':[], 'OL':[]}
    
    pct_ct = taxa_dict['total'][0]*pct_threshold
    pct_g = taxa_dict['total'][1]*pct_threshold
    pct_ol = taxa_dict['total'][2]*pct_threshold
    
    minimum_threshold = min([pct_ct, pct_g, pct_ol])
    
    taxa_list = list(all_taxa_dict.keys())
    taxa_list.sort(reverse=True)
    
    for taxa in taxa_list:
        abundances = all_taxa_dict[taxa]
        if taxa != 'total' and taxa != 'observed':
            s_ct, s_g, s_ol = sum(abundances[0]), sum(abundances[1]), sum(abundances[2]), 
            if max([s_ct, s_g, s_ol]) >= minimum_threshold:
                l_ct, l_g, l_ol = sum(return_log10(abundances[0], True)), sum(return_log10(abundances[1], True)), sum(return_log10(abundances[2], True))
                specific_taxa_dict['taxa'].append(taxa) 
                specific_taxa_dict['CT'].append(l_ct)
                specific_taxa_dict['G'].append(l_g)
                specific_taxa_dict['OL'].append(l_ol)

                outline = ('{}\t{}\t{}\t{}\n').format(taxa, l_ct, l_g, l_ol)
                logfile.write(outline)
    
    logfile.close()        
    outfile_name = ('taxonomic_abundance_{}_{}_{}_.pdf').format(prefix_name, pct_threshold, taxa_cutoff_name)
        
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        y=specific_taxa_dict['taxa'],
        x=specific_taxa_dict['CT'],
        marker=dict(color='rgba(0, 89, 255 0.5)', size=10),
        mode="markers",
        name="Sublittoral",
    ))
        
    fig.add_trace(go.Scatter(
        y=specific_taxa_dict['taxa'],
        x=specific_taxa_dict['G'],
        marker=dict(color='rgba(255, 165, 0, 0.5)', size=10),
        mode="markers",
        name="Intertidal",
    ))
    
    fig.add_trace(go.Scatter(
        y=specific_taxa_dict['taxa'],
        x=specific_taxa_dict['OL'],
        marker=dict(color='rgba(0, 0, 0, 0.5)', size=10),
        mode="markers",
        name="Supralittoral",
    ))

    fig.update_layout(title='Family Level Taxonomic Abundances',
                      xaxis_title="Log10 Taxa Abundance",
                      yaxis_title="Taxa",
                      font_size=10,
                      width=1500,
                      height=1500)
    
    fig.write_image(outfile_name)
        
for taxa_cutoff_name in rank_order:
    taxa_cutoff_num = convert_taxa_to_rank[taxa_cutoff_name]
    
    taxa_set = set()
    taxa_to_otu_dict = {}
    otu_to_taxa_dict = {}
    
    taxa_file = open(taxa_file_name)
    
    for line in taxa_file:
        if line[0]!='#':
            line = line.replace('"','')
            line = line.strip()
            otu = line.split('\t')[0]
            taxa = line.split('\t')[1]

            for each in prefixe:
                taxa = taxa.replace(each,'')

            taxa = taxa.replace(';','_')
            while taxa[-1] == '_':
                taxa = taxa[:-1]
                                   
            if taxa.count('_') >= taxa_cutoff_num:
                taxa_set.add(taxa)
                
                if taxa.count('_') > taxa_cutoff_num:
                    taxa_list = taxa.split('_')[:taxa_cutoff_num+1]
                    taxa = ''
                    for each in taxa_list:
                        taxa+=str(each)+'_'
                    
                    if taxa[-1] == '_':
                        taxa = taxa[:-1]
                                
                if taxa not in taxa_to_otu_dict:
                    taxa_to_otu_dict[taxa] = [otu]
                else:
                    taxa_to_otu_dict[taxa].append(otu)
                    
                if otu not in otu_to_taxa_dict:
                    otu_to_taxa_dict[otu] = taxa
    
                else:
                    print('err')
                
    taxa_file.close()
    
    otu_counts = build_otu_counts(feature_table_name)
    ct_cor_1, ct_cor_2, ct_cor_3, g_cor_1, g_cor_2, g_cor_3, ol_cor_1, ol_cor_2, ol_cor_3 = find_correction_value(otu_counts)
    
    taxa_to_counts = {}
    for taxa, otus in taxa_to_otu_dict.items():
        sub1, sub2, sub3 = 0, 0, 0
        int1, int2, int3 = 0, 0, 0
        sup1, sup2, sup3 = 0, 0, 0
        
        for otu in otus:
            if otu in otu_counts:
               osub1, osub2, osub3, oint1, oint2, oint3, osup1, osup2, osup3 = otu_counts[otu]
               sub1 += osub1
               sub2 += osub2
               sub3 += osub3
               
               int1 += oint1
               int2 += oint2
               int3 += oint3
               
               sup1 += osup1
               sup2 += osup2
               sup3 += osup3
               
        taxa_to_counts[taxa] = [sub1, sub2, sub3, int1, int2, int3, sup1, sup2, sup3]
    
    outfile = open(args.output_file+'_counts.tsv', 'w')
    header = ('#taxa\tP1\tP2\tP3\tS1\tS2\tS3\tV1\tV2\n')
    outfile.write(header)
    
    for taxa, cts in taxa_to_counts.items():
        sub1, sub2, sub3, int1, int2, int3, sup1, sup2, sup3 = cts
        outline = ('{taxa}\t{sub1}\t{sub2}\t{sub3}\t{int1}\t{int2}\t{int3}\t{sup1}\t{sup2}\t{sup3}\n').format(taxa=taxa, int1=int1, int2=int2, int3=int3, sub1=sub1, sub2=sub2, sub3=sub3, sup1=sup1, sup2=sup2, sup3=sup3)        
        outfile.write(outline)
    outfile.close()
        
    #store otu data
    taxa_dict = {'total':[0,0,0], 'observed':0}  
    
    taxa_raw_dict = {}
    
    for taxa, counts in taxa_to_counts.items():
        #
        sub1 = counts[0]
        sub2 = counts[1]
        sub3 = counts[2]
        raw_ct = [sub1, sub2, sub3]
        #
        int1 = counts[3]
        int2 = counts[4]  
        int3 = counts[5]
        raw_g = [int1, int2, int3]
        #
        sup1 = counts[6]
        sup2 = counts[7]  
        sup3 = counts[8]
        raw_ol = [sup1, sup2, sup3]
        #      
        if taxa not in taxa_raw_dict:
            taxa_raw_dict[taxa] = [raw_ct, raw_g, raw_ol]
        else:
            for index in range(3):
                taxa_raw_dict[taxa][0][index] += raw_ct[index]
                taxa_raw_dict[taxa][1][index] += raw_g[index]
                taxa_raw_dict[taxa][2][index] += raw_ol[index]
        #
        sub1 = counts[0]*ct_cor_1
        sub2 = counts[1]*ct_cor_2
        sub3 = counts[2]*ct_cor_3
        cor_ct = [sub1, sub2, sub3]
        #
        int1 = counts[3]*g_cor_1
        int2 = counts[4]*g_cor_2
        int3 = counts[5]*g_cor_3
        cor_g = [int1, int2, int3]
        #
        sup1 = counts[6]*ol_cor_1
        sup2 = counts[7]*ol_cor_2
        sup3 = counts[8]*ol_cor_3
        cor_ol = [sup1, sup2, sup3]
        #           
        if taxa not in taxa_dict:
            taxa_dict[taxa] = [cor_ct, cor_g, cor_ol]
        else:
            for index in range(3):
                taxa_dict[taxa][0][index] += cor_ct[index]
                taxa_dict[taxa][1][index] += cor_g[index]
                taxa_dict[taxa][2][index] += cor_ol[index]
                
        taxa_dict['total'][0] += sum(cor_ct)
        taxa_dict['total'][1] += sum(cor_g)
        taxa_dict['total'][2] += sum(cor_ol)
        
        taxa_dict['observed'] += len([s for s in cor_ct if s != 1]) + len([s for s in cor_g if s != 1]) + len([s for s in cor_ol if s != 1])
        
    header = ('taxa\tsub\tinter\tsup\tsub_bet\tint_bet\tsup_bet\tmin_pval\n')
    
    all_outfile_name = ('all_unnormalized_taxa_abundance_{}.tab').format(taxa_cutoff_name)
    all_outfile = open(all_outfile_name, 'w')
    all_outfile.write(header)
    
    for taxa, raw_taxa_array  in taxa_raw_dict.items():
        sub_s = sum(raw_taxa_array[0])
        int_s = sum(raw_taxa_array[1])
        sup_s = sum(raw_taxa_array[2])
        
        total_s = sub_s+int_s+sup_s
        pval = 1
        if total_s >= 100:
            if max((sub_s/total_s),(int_s/total_s),(sup_s/total_s)) >= 0.05:
                sub_stat = stats.binom_test(sub_s, n=total_s, p=0.333, alternative='two-sided')*len(taxa_raw_dict)
                int_stat = stats.binom_test(int_s, n=total_s, p=0.333, alternative='two-sided')*len(taxa_raw_dict)
                sup_stat = stats.binom_test(sup_s, n=total_s, p=0.333, alternative='two-sided')*len(taxa_raw_dict)
        
                if len([x for x in [sub_stat, int_stat, sup_stat] if x <= 0.05]) > 1:
                    pval = min([sub_stat, int_stat, sup_stat])
                
                
        outline = ('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(taxa, sub_s, int_s, sup_s, sub_stat, int_stat, sup_stat, pval)
        all_outfile.write(outline)

    all_outfile.close()
    
    all_outfile_name = ('all_normalized_taxa_abundance_{}.tab').format(taxa_cutoff_name)
    all_outfile = open(all_outfile_name, 'w')
    all_outfile.write(header)
    
    for taxa, taxa_array in taxa_dict.items():
        if taxa != 'total' and taxa != 'observed':
            sub_s = sum(taxa_array[0])
            int_s = sum(taxa_array[1])
            sup_s = sum(taxa_array[2])
            
            total_s = sub_s+int_s+sup_s
            pval = 1
            if total_s >= 100:
                if max((sub_s/total_s),(int_s/total_s),(sup_s/total_s))>=0.05:
                    sub_stat = stats.binom_test(sub_s, n=total_s, p=0.333, alternative='two-sided')*len(taxa_raw_dict)
                    int_stat = stats.binom_test(int_s, n=total_s, p=0.333, alternative='two-sided')*len(taxa_raw_dict)
                    sup_stat = stats.binom_test(sup_s, n=total_s, p=0.333, alternative='two-sided')*len(taxa_raw_dict)
            
                    if len([x for x in [sub_stat, int_stat, sup_stat] if x <= 0.05]) > 1:
                        pval = min([sub_stat, int_stat, sup_stat])
            
            outline = ('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(taxa, sub_s, int_s, sup_s, sub_stat, int_stat, sup_stat, pval)
            all_outfile.write(outline)
        
    all_outfile.close()
    
    #        
    taxa_outfile_name = ('site_specific_{}_enrichment.tab').format(taxa_cutoff_name)
    outfile = open(taxa_outfile_name, 'w')
    #taxa, uid, pval, p_1v2, p_1v3, p_1v4, p_2v3, p_2v4, p_3v4
    header = ('#taxa\tuid\tpval\tmedian_log_10_sub\tmedian_log_10_inter\tmedian_log_10_supra\tpval_sub_inter\tpval_sub_supra\tpval_inter_supra\tishow\n')
    outfile.write(header)
    
    plot_top10_taxa(taxa_dict, 'normalized', 0.01)
    plot_top10_taxa(taxa_raw_dict, 'un_normalized', 0.01)
    
    uid = 0
    
    max_obs = 0
    figure_dict = {}
    
    dotplot_dict = {}

    total_ct = taxa_dict['total'][0]    
    total_g = taxa_dict['total'][1]
    total_ol = taxa_dict['total'][2]
    
    total_observed = taxa_dict['observed']
    
    unique_dict = {}
        
    for taxa in taxa_set:                        
        #bonferroni_corrected_pvalue = 0.05
        if taxa in taxa_dict:
            #
            ct_r_array, ct_mean, ct_std, ct_set, ct_obs = return_deets(taxa_raw_dict[taxa][0], taxa_dict[taxa][0])
            g_r_array, g_mean, g_std, g_set, g_obs = return_deets(taxa_raw_dict[taxa][1], taxa_dict[taxa][1])
            ol_r_array, ol_mean, ol_std, ol_set, ol_obs = return_deets(taxa_raw_dict[taxa][2], taxa_dict[taxa][2])
            #
            ct_array, ct_mean, ct_std, ct_set, ct_obs = return_deets(taxa_dict[taxa][0], taxa_dict[taxa][0])
            g_array, g_mean, g_std, g_set, g_obs = return_deets(taxa_dict[taxa][1], taxa_dict[taxa][1])
            ol_array, ol_mean, ol_std, ol_set, ol_obs = return_deets(taxa_dict[taxa][2], taxa_dict[taxa][2])
            
            #unique function
            if ct_obs or g_obs or ol_obs:                
                if (sum(ct_r_array) >= 20 and sum(g_r_array) == 0) or (sum(ct_r_array) >= 20 and sum(ol_r_array) == 0):
                    ct = 0
                    for each in ct_r_array:
                        if each >= 10:
                            ct+=1
                            
                    if ct >= 2:
                        unique_dict[taxa]=[np.mean(return_log10(ct_array)), np.mean(return_log10(g_array)), np.mean(return_log10(ol_array))] 
                
                if (sum(g_r_array) >= 20 and sum(ct_r_array) == 0) or (sum(g_r_array) >= 20 and sum(ol_r_array) == 0): 
                    ct = 0
                    for each in g_r_array:
                        if each >= 10:
                            ct+=1
                            
                    if ct >= 2:
                        unique_dict[taxa]=[np.mean(return_log10(ct_array)), np.mean(return_log10(g_array)), np.mean(return_log10(ol_array))] 
                
                if (sum(ol_r_array) >= 20 and sum(g_r_array) == 0) or (sum(ol_r_array) >= 20 and sum(ct_r_array) == 0):
                    ct = 0
                    for each in ol_r_array:
                        if each >= 10:
                            ct+=1
                            
                    if ct >= 2:
                        unique_dict[taxa]=[np.mean(return_log10(ct_array)), np.mean(return_log10(g_array)), np.mean(return_log10(ol_array))] 
                
            #                                      
            if ct_obs or g_obs or ol_obs:
                p_1v2, p_1v3, p_2v3 = 0, 0, 0
                ct_sum, g_sum, ol_sum = sum(ct_array), sum(g_array), sum(ol_array)
                ct_less, g_less, ol_less = (total_ct - ct_sum), (total_g - g_sum), (total_ol - ol_sum)
                
                obs = np.array([[max(g_sum,1), g_less], [max(ct_sum,1), ct_less], [max(ol_sum,1), ol_less]])
                chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
                                              
                #stats_runmode = 'kruskal_wallis'
                
                if stats_runmode == 'kruskal_wallis':
                    pval = run_kruskal_3way(ct_set, g_set, ol_set)
                    pass_criteria = False
                    if pval <= 0.05:
                        pass_criteria = True
                      
                else:
                    p_1v2 = run_chi2(ct_sum, ct_less, g_sum, g_less)
                    p_1v3 = run_chi2(ct_sum, ct_less, ol_sum, ol_less)
                    p_2v3 = run_chi2(g_sum, g_less, ol_sum, ol_less)
                    pass_criteria = criteria(ct_array, g_array, ol_array, p_1v2, p_1v3, p_2v3, pval, taxa, pct_effect_size, pval_threshold)
                
                if pass_criteria:
                    log_ct_set = return_log10(ct_array)
                    log_g_set = return_log10(g_array)
                    log_ol_set = return_log10(ol_array)
                    
                    max_obs = test_max([log_ct_set, log_g_set, log_ol_set], max_obs)
                    
                    clog_ct_set = log_ct_set
                    clog_g_set = log_g_set
                    clog_ol_set = log_ol_set
                    
                    if (len(clog_ct_set)+len(clog_g_set)+len(clog_ol_set)) >= 3:
                        if clog_ol_set:
                            med_log_ol = np.median(clog_ol_set)
                        else:
                            med_log_ol = 0
                            
                        if clog_g_set:
                            med_log_g = np.median(clog_g_set)
                        else:
                            med_log_g = 0

                        if clog_ct_set:
                            med_log_ct = np.median(clog_ct_set)
                        else:
                            med_log_ct = 0
                        
                        if taxa not in figure_dict:
                            pass_dict['figure_dict']+=1
                            figure_dict[taxa] = [clog_ct_set, clog_g_set, clog_ol_set]
                        else:
                            print('error')
                            1/0
                        
                        is_how = simplify_enrichment(taxa, p_1v2, p_1v3, p_2v3, med_log_ct, med_log_g, med_log_ol)
                        
                        if is_how != 'complex':
                            outline = ('{taxa}\t{uid}\t{pval}\t{median_ct}\t{median_g}\t{median_ol}\t{pv12}\t{pv13}\t{pv23}\t{is_how}\n').format(taxa=taxa, uid=uid, pval=pval, median_ct=med_log_ct, median_g=med_log_g, median_ol=med_log_ol, pv12=p_1v2, pv13=p_1v3, pv23=p_2v3, is_how=is_how)
                            print(outline)
                            outfile.write(outline)
                            
                        uid += 1  
                        if taxa not in dotplot_dict:
                            dotplot_dict[taxa] = [med_log_ct, med_log_g, med_log_ol]
                            
                        else:
                            print('taxa conflict', taxa)
                           
    if unique_dict:
        import plotly.graph_objects as go
        
        taxa_by_site_list = []
        site_abundance_dict = {'CT':[], 'G':[], 'OL':[]}
        
        taxa_list = list(unique_dict.keys())
        taxa_list.sort(reverse=True)
            
        for taxa in taxa_list:
            sites = unique_dict[taxa]
            taxa_by_site_list.append(taxa)
            site_abundance_dict['CT'].append(sites[0])
            site_abundance_dict['G'].append(sites[1])
            site_abundance_dict['OL'].append(sites[2])
            
        outfile_name = ('_unique_{}_.pdf').format(taxa_cutoff_name)
            
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['CT'],
            marker=dict(color='rgba(0, 89, 255, 0.5)', size=10),
            mode="markers",
            name="Sub",
        ))
                
        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['G'],
            marker=dict(color='rgba(255, 165, 0, 0.5)', size=10),
            mode="markers",
            name="Inter",
        ))

        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['OL'],
            marker=dict(color='rgba(0, 0, 0, 0.5)', size=10),
            mode="markers",
            name="Sup",
        ))

        fig.update_layout(title='Site specific taxonomic differentials',
                          xaxis_title="Median Log10 Taxa Abundance",
                          yaxis_title="Taxa",
                          font_size=10,
                          width=1500,
                          height=1500)
        
        fig.write_image(outfile_name)     
        
        
    if figure_dict:
        #build_x_y()          
        for taxa, set_list in figure_dict.items():
            
            global_x_data = []
            global_y_data = []
            global_colors = []
            
            if '/' in taxa:
                taxa = taxa.replace('/','_or_')
            
            ct_tag = ('{}_{}_Sub').format(taxa_cutoff_name, taxa)
            g_tag = ('{}_{}_Inter').format(taxa_cutoff_name, taxa)
            ol_tag = ('{}_{}_Sup').format(taxa_cutoff_name, taxa)
             
            x_data = ct_tag, g_tag, ol_tag
            
            ct_set = mod_null_set(set_list[0])            
            g_set = mod_null_set(set_list[1])
            ol_set = mod_null_set(set_list[2])
                    
            outfile.close()
                                                           
            y_data = ct_set, g_set, ol_set
            
            colors = 'rgba(0, 89, 255 0.5)', 'rgba(255, 165, 0, 0.5)', 'rgba(0, 0, 0, 0.5)',
            
            global_x_data += x_data
            global_y_data += y_data
            global_colors += colors
            
            fig = go.Figure()
            
            outfile_name = ('site_specific_{}_{}_enrichment.pdf').format(taxa_cutoff_name, taxa)
            print(outfile_name)
            
            for xd, yd, cls in zip(global_x_data, global_y_data, global_colors):
                    fig.add_trace(go.Box(
                        y=yd,
                        name=xd,
                        boxpoints='all',
                        notched=True,
                        jitter=0.5,
                        whiskerwidth=0.2,
                        fillcolor=cls,
                        line_color=cls,
                        marker_size=5,
                        line_width=1,
                        showlegend=False)
                    )
                    
            fig.update_layout(
                title=taxa,
                xaxis_title="Sample Site",
                yaxis_title="Log10(Relative Taxa Abundance)",

            )

            fig.write_image(outfile_name)
    
    if dotplot_dict:
        import plotly.graph_objects as go
        
        taxa_by_site_list = []
        site_abundance_dict = {'CT':[], 'G':[], 'OL':[]}
        
        taxa_list = list(dotplot_dict.keys())
        taxa_list.sort(reverse=True)
            
        for taxa in taxa_list:
            sites = dotplot_dict[taxa]
            taxa_by_site_list.append(taxa)
            site_abundance_dict['CT'].append(sites[0])
            site_abundance_dict['G'].append(sites[1])
            site_abundance_dict['OL'].append(sites[2])
            
        outfile_name = ('_site_specific_spread_{}_.pdf').format(taxa_cutoff_name)
            
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['CT'],
            marker=dict(color='rgba(0, 89, 255 0.5)', size=10),
            mode="markers",
            name="Sublittoral",
        ))
                
        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['G'],
            marker=dict(color='rgba(255, 165, 0, 0.5)', size=10),
            mode="markers",
            name="Intertidal",
        ))

        fig.add_trace(go.Scatter(
            y=taxa_by_site_list,
            x=site_abundance_dict['OL'],
            marker=dict(color='rgba(0, 0, 0, 0.5)', size=10),
            mode="markers",
            name="Supralittoral",
        ))

        fig.update_layout(title='Site specific taxonomic differentials',
                          xaxis_title="Median Log10 Taxa Abundance",
                          yaxis_title="Taxa",
                          font_size=10,
                          width=1500,
                          height=1500)
        
        fig.write_image(outfile_name)     