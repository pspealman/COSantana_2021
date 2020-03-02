# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 09:29:45 2019

Novembro - for identifying taxa enrichments from Qiime generated data

python novembro.py -f feature-table.biom.txt -t taxonomy.tsv -s silva -o taxa_counts.tab

NB: feature-table.biom.txt should be generated from converting the qiime feature_table.biom file to a tsv
#   biom convert -i feature-table.biom -o feature-table.biom.txt --to-tsv

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
args = parser.parse_args()

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
    
otu_file = open(feature_table_name)

otu_counts = {}

for line in otu_file:
    if line[0]!='#' and 'Feature ID' not in line:
        line = line.strip()
        otu = line.split('\t')[0].strip()
#        
        submerged_1 = float(line.split('\t')[1])
        submerged_2 = float(line.split('\t')[2])       
        submerged_3 = float(line.split('\t')[3])
        #
        intertidal_1 = float(line.split('\t')[4])
        intertidal_2 = float(line.split('\t')[5])  
        intertidal_3 = float(line.split('\t')[6])
        #
        seco_1 = float(line.split('\t')[7])
        seco_2 = float(line.split('\t')[8])
        seco_3 = float(line.split('\t')[9])
#        #
        otu_counts[otu] = [submerged_1, submerged_2, submerged_3, intertidal_1, intertidal_2, intertidal_3, seco_1, seco_2, seco_3]
#        
otu_file.close()

taxa_to_otu_dict = {}

taxa_file = open(taxa_file_name)

for line in taxa_file:
    if line[0]!='#' and 'Feature ID' not in line:
        line = line.replace('"','')
        line = line.strip()
        otu = line.split('\t')[0]
        taxa = line.split('\t')[1]
        taxa = taxa.replace(' ', '_')
        for each in prefixe:
            taxa = taxa.replace(each,'')

        taxa = taxa.replace(';','_')
        while taxa[-1] == '_':
            taxa = taxa[:-1]
                                                
        if taxa not in taxa_to_otu_dict:
            otu_set = set()
            otu_set.add(otu)
            taxa_to_otu_dict[taxa] = otu_set
        else:
            taxa_to_otu_dict[taxa].add(otu)
                
taxa_file.close()

taxa_to_counts = {}

for taxa, otus in taxa_to_otu_dict.items():
    submerged_1, submerged_2, submerged_3 = 0, 0, 0
    intertidal_1, intertidal_2, intertidal_3 = 0, 0, 0
    seco_1, seco_2, seco_3 = 0, 0, 0 
    
    for otu in otus:
       nsubmerged_1, nsubmerged_2, nsubmerged_3, nintertidal_1, nintertidal_2, nintertidal_3, nseco_1, nseco_2, nseco_3 = otu_counts[otu]
       submerged_1 += nsubmerged_1
       submerged_2 += nsubmerged_2
       submerged_3 += nsubmerged_3
       intertidal_1 += nintertidal_1
       intertidal_2 += nintertidal_2
       intertidal_3 += nintertidal_3
       seco_1 += nseco_1
       seco_2 += nseco_2
       seco_3 += nseco_3
       
    taxa_to_counts[taxa] = [submerged_1, submerged_2, submerged_3, intertidal_1, intertidal_2, intertidal_3, seco_1, seco_2, seco_3]

outfile = open(args.output_file, 'w')
header = ('#taxa\tP1_1\tP1_2\tP1_3\tP2_1\tP2_2\tP2_3\tP3_1\tP3_2\tP3_3\n')
outfile.write(header)

for taxa, cts in taxa_to_counts.items():
    submerged_1, submerged_2, submerged_3, intertidal_1, intertidal_2, intertidal_3, seco_1, seco_2, seco_3 = cts
    outline = ('{taxa}\t{p11}\t{p12}\t{p13}\t{p21}\t{p22}\t{p23}\t{p31}\t{p32}\t{p33}\n').format(taxa=taxa, p11=submerged_1, p12=submerged_2, p13=submerged_3, p21=intertidal_1, p22=intertidal_2, p23=intertidal_3, p31=seco_1, p32=seco_2, p33=seco_3)        
    outfile.write(outline)
outfile.close()

#
rank_order = ['species','genus','family','order','class','phylum','kingdom']
convert_taxa_to_rank = {'kingdom':0, 'phylum':1, 'class':2, 'order':3, 'family':4, 'genus':5, 'species': 6}

filter_previous = []
pass_dict = {'criteria':0, 'pval':0, 'p_pval':0, 'pass_set':0, 'figure_dict':0}

for taxa_cutoff_name in rank_order:
    taxa_cutoff_num = convert_taxa_to_rank[taxa_cutoff_name]
    
    taxa_set = set()
    taxa_to_otu_dict = {}
    otu_to_taxa_dict = {}
    
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
                
            if taxa not in filter_previous:                    
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
    
    otu_file = open(feature_table_name)
    
    seco_1, seco_2, seco_3 = 0, 0, 0 
    intertidal_1, intertidal_2, intertidal_3 = 0, 0, 0
    submerged_1, submerged_2, submerged_3 = 0, 0, 0
    
    for line in otu_file:
        if line[0]!='#':
            submerged_1 += float(line.split('\t')[1])
            submerged_2 += float(line.split('\t')[2])       
            submerged_3 += float(line.split('\t')[3])
            #
            intertidal_1 += float(line.split('\t')[4])
            intertidal_2 += float(line.split('\t')[5])  
            intertidal_3 += float(line.split('\t')[6])
            #
            seco_1 += float(line.split('\t')[7])
            seco_2 += float(line.split('\t')[8])
            seco_3 += float(line.split('\t')[9])
            #
            
    otu_file.close()
        
    global_min = min(sum([seco_1,seco_2,seco_3]),
                     sum([intertidal_1,intertidal_2,intertidal_3]),
                     sum([submerged_1,submerged_2,submerged_3]))
    
    #Define number of observations per site for the purposes of downsampling
    seco_cor_1, seco_cor_2, seco_cor_3 = global_min/seco_1, global_min/seco_2, global_min/seco_3 
    intertidal_cor_1, intertidal_cor_2, intertidal_cor_3 = global_min/intertidal_1, global_min/intertidal_2, global_min/intertidal_3
    submerged_cor_1, submerged_cor_2, submerged_cor_3 = global_min/submerged_1, global_min/submerged_2, global_min/submerged_3  
    
    def obs_counter(list_obs_ct):
        obs = 0
        for obs_ct in list_obs_ct:
            if obs_ct > 0:
                obs+=1
        
        if obs > 1:
            return(True)
        else:
            return(False)
        
    #store otu data
    taxa_dict = {'total':[0,0,0], 'observed':0}  
    
    taxa_raw_dict = {}
    
    otu_file = open(feature_table_name)
    
    for line in otu_file:
        if line[0]!='#':
            otu = line.split('\t')[0]
            
            if otu in otu_to_taxa_dict:
                taxa = otu_to_taxa_dict[otu]
                                
                submerged_1 = float(line.split('\t')[1])
                submerged_2 = float(line.split('\t')[2])
                submerged_3 = float(line.split('\t')[3])
                raw_submerged = [submerged_1, submerged_2, submerged_3]
                 
                intertidal_1 = float(line.split('\t')[4])
                intertidal_2 = float(line.split('\t')[5])
                intertidal_3 = float(line.split('\t')[6])
                raw_intertidal = [intertidal_1, intertidal_2, intertidal_3]
                                #
                seco_1 = float(line.split('\t')[7])
                seco_2 = float(line.split('\t')[8])
                seco_3 = float(line.split('\t')[9])
                raw_seco = [seco_1, seco_2, seco_3]
  
                if taxa not in taxa_raw_dict:
                    taxa_raw_dict[taxa] = [raw_seco, raw_intertidal, raw_submerged]
                else:
                    taxa_raw_dict[taxa][0] += raw_seco
                    taxa_raw_dict[taxa][1] += raw_intertidal
                    taxa_raw_dict[taxa][2] += raw_submerged

                submerged_1 = (float(line.split('\t')[1])*submerged_cor_1)
                submerged_2 = (float(line.split('\t')[2])*submerged_cor_2)
                submerged_3 = (float(line.split('\t')[3])*submerged_cor_3)
                submerged = [submerged_1, submerged_2, submerged_3]
                #
                intertidal_1 = (float(line.split('\t')[4])*intertidal_cor_1)
                intertidal_2 = (float(line.split('\t')[5])*intertidal_cor_2)
                intertidal_3 = (float(line.split('\t')[6])*intertidal_cor_3)
                intertidal = [intertidal_1, intertidal_2, intertidal_3]
                #
                seco_1 = (float(line.split('\t')[7])*seco_cor_1)
                seco_2 = (float(line.split('\t')[8])*seco_cor_2)
                seco_3 = (float(line.split('\t')[9])*seco_cor_3)
                seco = [seco_1, seco_2, seco_3]
                #           
                if taxa not in taxa_dict:
                    taxa_dict[taxa] = [seco, intertidal, submerged]
                else:
                    taxa_dict[taxa][0] += seco
                    taxa_dict[taxa][1] += intertidal
                    taxa_dict[taxa][2] += submerged
                    
                taxa_dict['total'][0] += sum(seco)
                taxa_dict['total'][1] += sum(intertidal)
                taxa_dict['total'][2] += sum(submerged)
                
                taxa_dict['observed'] += len([s for s in seco if s != 1]) + len([s for s in intertidal if s != 1]) + len([s for s in submerged if s != 1])
    otu_file.close()

    
    all_outfile_name = ('all_taxa_abundance_{}.tab').format(taxa_cutoff_name)
    all_outfile = open(all_outfile_name, 'w')
    for taxa, raw_taxa_array  in taxa_raw_dict.items():
        print(raw_taxa_array)
        1/0
        outline = ('{}\t{}\t{}\t{}\n').format(taxa, sum(raw_taxa_array[0]), sum(raw_taxa_array[1]), sum(raw_taxa_array[2]))
        all_outfile.write(outline)
    all_outfile.close()

    #        
    taxa_outfile_name = ('site_specific_{}_enrichment.tab').format(taxa_cutoff_name)
    outfile = open(taxa_outfile_name, 'w')
    header = ('#taxa\tuid\tchi_pval\tMWU_pval_seco_intertidal\tMWU_pval_intertidal_submerged\tMWU_pval_submerged_seco\n')
    outfile.write(header)
        
    uid = 0
    
    max_obs = 0
    figure_dict = {}
    
    def criteria(seco_set, intertidal_set, submerged_set, p_sc, p_is, p_ss, pval, taxa, pct_effect_size=0.05):
        pass_set = []
        log_fold_diff = []
        
        global pass_dict
        pass_dict['criteria']+=1
        
        if pval <= 0.05:
            pass_dict['pval']+=1
            max_effect_size = (1000+pct_effect_size*1000)/1000
            min_effect_size = (1000-pct_effect_size*1000)/1000
            
            seco_mean = np.mean(seco_set)
            intertidal_mean = np.mean(intertidal_set)
            submerged_mean = np.mean(submerged_set)
                    
            if seco_mean == 0:
                seco_mean = 1
            if intertidal_mean == 0:
                intertidal_mean = 1
            if submerged_mean == 0:
                submerged_mean = 1
                                
            w_sc = (seco_mean/intertidal_mean)
            w_is = (intertidal_mean/submerged_mean)
            w_ss = (submerged_mean/seco_mean)
            
            log_fold_diff = [w_sc, w_is, w_ss]
                     
            if (p_sc <= 0.05) or (p_is <= 0.05) or (p_ss <= 0.05):
                pass_dict['p_pval'] += 1
            
            if (p_sc <= 0.05) and (w_sc >= max_effect_size or w_sc <= min_effect_size):
                pass_set.append('sc')
                print(w_sc)
                    
            if (p_is <= 0.05) and (w_is >= max_effect_size or w_is <= min_effect_size):
                pass_set.append('is')
                print(w_is)
            
            if (p_ss <= 0.05) and (w_ss >= max_effect_size or w_ss <= min_effect_size):
                pass_set.append('ss')
                print(w_ss)
                
        if len(pass_set) > 0:
            pass_dict['pass_set']+=1
            return(True, log_fold_diff, pass_set)
        else:
            return(False, log_fold_diff, pass_set)
    
    def test_max(set_list, max_obs):
        for each_set in set_list:
            if max(each_set) >= max_obs:
                max_obs = max(each_set)
        return(max_obs)
    
    def return_log10(each_set):
        new_set = []
        for each_obs in each_set:
            if each_obs <= 0:
                each_obs = 1
            else:
                each_obs = np.log10(each_obs)

            new_set.append(each_obs)
        
        return(new_set)
        
    figure_dict = {}
    
    total_seco = taxa_dict['total'][0]
    total_intertidal = taxa_dict['total'][1]
    total_submerged = taxa_dict['total'][2]
    
    total_observed = taxa_dict['observed']
        
    for taxa in taxa_set:
        if taxa in taxa_dict and taxa in taxa_raw_dict:
            seco_array, intertidal_array, submerged_array = taxa_dict[taxa]
            #
            seco_mean = np.mean(seco_array)
            seco_std = np.std(seco_array)
            intertidal_mean = np.mean(intertidal_array)
            intertidal_std = np.std(intertidal_array)       
            submerged_mean = np.mean(submerged_array)
            submerged_std = np.std(submerged_array)
            #
            raw_taxa_array = taxa_raw_dict[taxa]
            seco_set = raw_taxa_array[0]
            intertidal_set = raw_taxa_array[1]
            submerged_set = raw_taxa_array[2]
                                        
            seco_obs = obs_counter(seco_set)
            intertidal_obs = obs_counter(intertidal_set)
            submerged_obs = obs_counter(submerged_set)
                            
            if seco_obs or intertidal_obs or submerged_obs:
                p_sc = 0
                p_ss = 0
                p_is = 0
                
                obs = np.array([[max(seco_mean,1), (total_seco - sum(seco_array))], [max(intertidal_mean,1), (total_intertidal - sum(intertidal_array))], [max(submerged_mean,1), (total_submerged - sum(submerged_array))]])
            
                chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
                                
                if sum(seco_set) > 30 or sum(intertidal_set) > 30:
                    _w, p_sc = stats.mannwhitneyu(seco_set, intertidal_set)
                        
                if sum(seco_set) > 30 or sum(submerged_set) > 30:
                    _w, p_ss = stats.mannwhitneyu(seco_set, submerged_set)
        
                if sum(intertidal_set) > 30 or sum(submerged_set) > 30:
                    _w, p_is = stats.mannwhitneyu(intertidal_set, submerged_set)
                
                pass_criteria, log_fold_difference, pass_set = criteria(seco_array, intertidal_array, submerged_array, p_sc, p_is, p_ss, pval, taxa)
        
                if pass_criteria:
                    log_seco_set = return_log10(seco_array)
                    log_intertidal_set = return_log10(intertidal_array)
                    log_submerged_set = return_log10(submerged_array)
                    
                    max_obs = test_max([log_seco_set,log_intertidal_set,log_submerged_set], max_obs)
                    
                    clog_seco_set = [s for s in log_seco_set if s != 1]
                    clog_intertidal_set = [s for s in log_intertidal_set if s != 1]
                    clog_submerged_set = [s for s in log_submerged_set if s != 1]
                    
                    if (len(clog_seco_set)+len(clog_intertidal_set)+len(clog_submerged_set)) >= (10):
                        if taxa not in figure_dict:
                            pass_dict['figure_dict']+=1
                            figure_dict[taxa] = [clog_seco_set, clog_intertidal_set, clog_submerged_set]
                        else:
                            print('error')
                            1/0
            
                        outline = ('{}\t{}\t{}\t{}\t{}\t{}\n').format(taxa, uid, pval, p_sc, p_is, p_ss)
                        outfile.write(outline)
                        
                        uid += 1  
        
    if figure_dict:
        for taxa, set_list in figure_dict.items():
            
            filter_previous.append(taxa)
            
            if '/' in taxa:
                taxa = taxa.replace('/','_or_')
            
            seco_tag = ('{}_{}_Seco').format(taxa_cutoff_name, taxa)
            intertidal_tag = ('{}_{}_Intertidal').format(taxa_cutoff_name, taxa)
            submerged_tag = ('{}_{}_Submerged').format(taxa_cutoff_name, taxa)
            x_data = seco_tag, intertidal_tag,submerged_tag
            
            seco_set = set_list[0]
            intertidal_set = set_list[1]
            submerged_set = set_list[2]
                    
            outfile.close()
            
            global_x_data = []
            global_y_data = []
            global_colors = []
                                               
            y_data = seco_set, intertidal_set, submerged_set
            
            colors = 'rgba(255, 144, 14, 0.5)', 'rgba(44, 160, 101, 0.5)', 'rgba(93, 164, 214, 0.5)'
            
            global_x_data += x_data
            global_y_data += y_data
            global_colors += colors
            
        fig = go.Figure()
        
        outfile_name = ('site_specific_{}_enrichment.pdf').format(taxa_cutoff_name)
        
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
                
        fig.show()
        fig.write_image(outfile_name)
