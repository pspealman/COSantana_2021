# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 23:03:48 2020

sigilo - for the generation and analysis of PICRUSt2 output

python sigilo --generate_heatmap -i pred_metagenome_unstrat.tsv -sig significant_objects_file -o metagenome

python sigilo --find_enrichment -c pred_metagenome_contrib.tsv -t taxonomy.tsv -p metabolism_KOs/ -select select_taxa.tsv -pct 0.1 -pval 0.05 -o metagenome_contrib

@author: pspea
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-hm',"--generate_heatmap", action='store_true')
parser.add_argument('-i',"--input_abundance_file")
parser.add_argument('-sig',"--significant_objects_file")

parser.add_argument('-en',"--find_enrichment", action='store_true')
parser.add_argument('-c',"--contrib_file")
parser.add_argument('-t',"--taxonomy")
parser.add_argument('-select', '--select_taxa_list')
parser.add_argument('-pct', '--pct_threshold')
parser.add_argument('-pval', '--pval_threshold')
parser.add_argument('-p',"--path_to_ko_files")

parser.add_argument('-o',"--output_file")

args = parser.parse_args()

#import packages:
import numpy as np
import plotly.graph_objects as go
import scipy.stats as stats

def load_sig_obj():
    so_file = args.significant_objects
    sig_object_set = set()
    
    for line in so_file:
        if (line.split('\t')[0]):
            sig_object_set.add(line.split('\t')[0])
    so_file.close()
    
    return(sig_object_set)

def parse_line(line, runmode='log'):
    if runmode == 'log':
        KO = line.split('\t')[0].strip()
        #
        P1_1 = np.log(max(float(line.split('\t')[1]),1))
        P1_2 = np.log(max(float(line.split('\t')[2]),1))
        P1_3 = np.log(max(float(line.split('\t')[3]),1))
        #
        P2_1 = np.log(max(float(line.split('\t')[4]),1))
        P2_2 = np.log(max(float(line.split('\t')[5]),1))
        P2_3 = np.log(max(float(line.split('\t')[6]),1))        
        #
        P3_1 = np.log(max(float(line.split('\t')[7]),1))
        P3_2 = np.log(max(float(line.split('\t')[8]),1))
        P3_3 = np.log(max(float(line.split('\t')[9]),1))
        #
        description = line.split('\t')[0]
    else:
        KO = line.split('\t')[0].strip()
        #
        P1_1 = float(line.split('\t')[1])
        P1_2 = float(line.split('\t')[2])
        P1_3 = float(line.split('\t')[3])
        #
        P2_1 = float(line.split('\t')[4])
        P2_2 = float(line.split('\t')[5])
        P2_3 = float(line.split('\t')[6])        
        #
        P3_1 = float(line.split('\t')[7])
        P3_2 = float(line.split('\t')[8])
        P3_3 = float(line.split('\t')[9])
        #
        description = line.split('\t')[0]
    
    return(KO, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, P3_1, P3_2, P3_3, description)

def map_kos(path):
    import glob
    ko_to_object_dict = {}
    metabolism_dict = {}
    # examples of ordered lists from Alexi
    #ordered_list = ["K00886","K01007","K01715","K01849","K01959","K09709","K14466","K15016","K17067","K18472","K00131","K05774","K07404","K00852","K01625","K00874","K01455","K00370","K00371","K00374","K00956","K00958"]
    #ordered_list = ["K01055","K18249","K18248"]
    ordered_list = []
    
    ko_list = [f for f in glob.glob(path + "**/*.txt", recursive=True)]
       
    for each_file in ko_list:
        #print(each_file)
        infile = open(each_file)
        
        for line in infile:
            if line[0]=='#':
                name = line.split('=')[1].strip()
                metabolism_dict[name] = []

            else:
                line = line.strip()

                if '\t' in line:
                    metabolism_dict[name].append(line.split('\t')[0])
                    ordered_list.append(line.split('\t')[0])
                
                    ko, function = line.split('\t')
                
                    if ko not in ko_to_object_dict:
                        ko_to_object_dict[ko] = {function}
                    else:
                        ko_to_object_dict[ko].add(function)
                
        infile.close()
        
    return(metabolism_dict, ordered_list, ko_to_object_dict)
   
def go_master_heatmap(ko_master_dict):
    global max_value
    global ko_dict
    
    heatmap_name = ('{}_master_heatmap.pdf').format(args.output_file)
    
    title_is = ('Heatmap of Significant KOs').format()

    site_list = ["Sub_1", "Sub_2", "Sub_3", "Int_1", "Int_2", "Int_3", "Sec_1", "Sec_2", "Sec_3"]
    ko_list = []
    value_array = []
    
    for ko, data in ko_master_dict.items():
        #data=ko_master_dict[ko]
        print(ko, data)
        value_array.append(data)
        ko_line = ('{}: {}').format(ko, ko_dict[ko])
        ko_list.append(ko_line)

    fig = go.Figure(data=go.Heatmap(
                    z=value_array,
                    x=site_list,
                    y=ko_list,
                    hoverongaps = False,
                    zmin=0, zmax=max_value))
    
    fig.update_layout(
            title=title_is,
            autosize=False,
            width=3600,
            height=3600,
            margin=go.layout.Margin(
                    l=50,
                    r=50,
                    b=100,
                    t=100,
                    pad=4
                    )
            )

    fig.show()
    fig.write_image(heatmap_name)
        
def go_median_heatmap(ko_master_dict):
    global max_value
    global ko_dict
    
    heatmap_name = ('{}_median_heatmap.pdf').format(args.output_file)
    
    title_is = ('Heatmap of Significant KOs').format()

    site_list = ["Sub", "Intertidal", "Supra"]
    ko_list = []
    value_array = []
    
    for ko, data in ko_master_dict.items():
        data=ko_master_dict[ko]
        print(ko, data)
        value_array.append(data)
        ko_line = ('{}: {}').format(ko, ko_dict[ko])
        ko_list.append(ko_line)

    fig = go.Figure(data=go.Heatmap(
                    z=value_array,
                    x=site_list,
                    y=ko_list,
                    hoverongaps = False,
                    zmin=0, zmax=max_value))
    fig.update_layout(
            title=title_is,
            autosize=False,
            width=3600,
            height=3600,
            margin=go.layout.Margin(
                    l=50,
                    r=50,
                    b=100,
                    t=100,
                    pad=4
                    )
            )
    
    fig.show()
    fig.write_image(heatmap_name)
    
def with_sig_object():
        
    sig_object_set = load_sig_obj()
    ko_set = set()
    
    for pwy in sig_object_set:
        ko_set.add(pwy)
        
    abun_val_file = open(args.input_abundance_file)
    
    ko_dict = {}
    master_ko_dict = {}
    median_ko_dict = {}
    max_value = 0
    
    for line in abun_val_file:
        #KO	P1.1	P1.2	P1.3	P2.1	P2.2	P2.3	P3.1	P3.2	P3.3	description
        if line[0]!='#':
            KO, p1_1, p1_2, p1_3, p2_1, p2_2, p2_3, p3_1, p3_2, p3_3, description = parse_line(line)
            if KO in ko_set:
                max_value  = max([p1_1, p1_2, p1_3, p2_1, p2_2, p2_3, p3_1, p3_2, p3_3, max_value])
                        
                if KO not in ko_dict:
                    ko_dict[KO]=description.strip()
                    master_ko_dict[KO]=[p1_1, p1_2, p1_3, p2_1, p2_2, p2_3, p3_1, p3_2, p3_3]
                    median_ko_dict[KO]=[np.mean([p1_1, p1_2, p1_3]),np.mean([p2_1, p2_2, p2_3]),np.mean([p3_1, p3_2, p3_3])]
                else:
                    print('Error: Duplicate KO Identified')
                    quit()
    
    abun_val_file.close()
    
    go_master_heatmap(master_ko_dict)
    
    go_median_heatmap(median_ko_dict)


def without_sig_object():
    abun_val_file = open(args.input_abundance_file)
    
    ko_dict = {}
    master_ko_dict = {}
    median_ko_dict = {}

    max_value = 0
    
    for line in abun_val_file:
        #KO	P1.1	P1.2	P1.3	P2.1	P2.2	P2.3	P3.1	P3.2	P3.3	description
        if line[0]!='#':
            KO, p1_1, p1_2, p1_3, p2_1, p2_2, p2_3, p3_1, p3_2, p3_3, description = parse_line(line, 'nolog')
            
            obs = np.array([[p1_1, p2_1, p3_1], [p1_2, p2_2, p3_2], [p1_3, p2_3, p3_3]])
            chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
            
            if pval <= 0.05:
                KO, p1_1, p1_2, p1_3, p2_1, p2_2, p2_3, p3_1, p3_2, p3_3, description = parse_line(line)
                max_value  = max([p1_1, p1_2, p1_3, p2_1, p2_2, p2_3, p3_1, p3_2, p3_3, max_value])
                
                if KO not in ko_dict:
                    ko_dict[KO]=description.strip()
                    
                    master_ko_dict[KO]=[p1_1, p1_2, p1_3, p2_1, p2_2, p2_3, p3_1, p3_2, p3_3]
                    median_ko_dict[KO]=[np.mean([p1_1, p1_2, p1_3]),np.mean([p2_1, p2_2, p2_3]),np.mean([p3_1, p3_2, p3_3])]
                    print(KO, np.mean([p1_1, p1_2, p1_3]),np.mean([p2_1, p2_2, p2_3]),np.mean([p3_1, p3_2, p3_3]))
                else:
                    print('Error: Duplicate KO Identified')
                    quit()
    
    abun_val_file.close()
    
    go_master_heatmap(master_ko_dict)
    
    go_median_heatmap(median_ko_dict)
    

### are single taxa enriched in KOs
def parse_taxonomy(taxonomy_name):
    taxonomy = open(taxonomy_name)
    
    taxa_dict = {}
    
    for line in taxonomy:
        t_id = line.split('\t')[0]
        taxa = line.split('\t')[1]
                        
        if ';D_7__' in taxa:
            taxa = taxa.split(';D_7__')[0]
                
        taxa_dict[t_id]=taxa
        
    taxonomy.close()
    
    return(taxa_dict)

def parse_contrib(contrib_name):
    contrib =open(contrib_name)
    ko_dict = {}
    taxonomic_level_total = {}
    ko_tfa_dict = {}
    
    for line in contrib:
        ko = line.split('\t')[1]
        
        if line.split('\t')[1] in ordered_list:
            tid = line.split('\t')[2]
            
            taxon_function_abun = float(line.split('\t')[6])

            taxon = taxa_dict[tid]
                            
            while taxon.count(';') > 0:
                level = taxon.count(';')
                
                if level not in taxonomic_level_total:
                    taxonomic_level_total[level] = {ko:[taxon_function_abun]}
                else:
                    if ko not in taxonomic_level_total[level]:
                        taxonomic_level_total[level][ko] = [taxon_function_abun]
                    else:
                        taxonomic_level_total[level][ko].append(taxon_function_abun)
              
                if ko not in ko_dict:
                    ko_dict[ko] = {taxon:taxon_function_abun}
                else:
                    if taxon not in ko_dict[ko]:
                        ko_dict[ko][taxon] = taxon_function_abun
                    else:
                        ko_dict[ko][taxon] += taxon_function_abun
                    
                if ko not in ko_tfa_dict:
                    ko_tfa_dict[ko] = [taxon_function_abun]
                else:
                    ko_tfa_dict[ko].append(taxon_function_abun)
                
                taxon = taxon.rsplit(';',1)[0]
                     
    contrib.close()
    
    return(ko_dict, taxonomic_level_total, ko_tfa_dict)

def check_metabolism(ko):
    global metabolism_dict
    meta_string = ''
    for meta, ko_dict in metabolism_dict.items():
        if ko in ko_dict:
            meta_string += meta + ','
    return(meta_string)
    
def apply_first_round(pct_threshold=0.2, pval_threshold=0.05):
    global ko_dict, taxonomic_level_total, ko_tfa_dict, taxa_gfc
    
    correction = len(ko_dict)

    second_round = {}

    for ko, _taxon_function_abun in ko_tfa_dict.items():
        if ko in ko_dict:
            #taxons = 
            for taxon, taxa_fa in ko_dict[ko].items():
                
                level = taxon.count(';')
                taxon_function_abun = sum(taxonomic_level_total[level][ko])
                num_taxon_function_abun = len(taxonomic_level_total[level][ko])
                                   
                if (taxa_fa/taxon_function_abun) >= pct_threshold:
                    meta_string = check_metabolism(ko)
                    gfc = 1 #taxa_gfc[taxon] 
                    
                    pval = stats.binom_test(taxa_fa, n=(taxon_function_abun), p=(gfc/num_taxon_function_abun))*correction
                    
                    if pval <= pval_threshold:
                        meta_string = check_metabolism(ko)
                        
                        if taxon not in second_round:
                            second_round[taxon] = {}
                            
                        if meta_string not in second_round[taxon]:
                            second_round[taxon][meta_string] = set()

                        second_round[taxon][meta_string].add(ko)
                            
    return(second_round)


def apply_second_round(measures=3, pct_threshold=0.2, pval_threshold=0.05, select_taxa=False):
    global ko_dict, taxonomic_level_total, ko_tfa_dict, taxa_gfc, second_round, ko_to_object_dict
    
    correction = len(ko_dict)
    
    if select_taxa:
        outname = ('{}.tsv').format(select_taxa.rsplit('__',1)[1])
    else:
        outname = ('global_criteria_{}_{}.tsv').format(pct_threshold, pval_threshold)
    out_temp = open(outname,'w')
       
    for ko, _taxon_function_abun in ko_tfa_dict.items():
        if ko in ko_dict:
            for taxon, taxa_fa in ko_dict[ko].items():
                process = True
                
                if select_taxa:
                    process = False
                    if select_taxa in taxon:
                        process = True
                
                if process:
                    meta_string = check_metabolism(ko)
                    if taxon in second_round:
                        if meta_string in second_round[taxon]:
                            if len(second_round[taxon][meta_string]) >= measures:
                                level = taxon.count(';')
                                gfc = 1 #taxa_gfc[taxon]
                                
                                taxon_function_abun = sum(taxonomic_level_total[level][ko])
                                num_taxon_function_abun = len(taxonomic_level_total[level][ko])
                                
                                if (taxa_fa/taxon_function_abun) >= pct_threshold:
                                    pval = stats.binom_test(taxa_fa, n=(taxon_function_abun), p=(gfc/num_taxon_function_abun))*correction
                                    
                                    if pval <= 0.05:
                                        if (taxa_fa/taxon_function_abun) >= pval_threshold:                                    
                                            meta_string = check_metabolism(ko)
                                            function = ko_to_object_dict[ko]
                                            outline = ('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(meta_string, level, taxon, ko, function, (taxa_fa/taxon_function_abun), pval, taxa_fa, taxon_function_abun)
                                            print(outline)
                                            out_temp.write(outline)                                           
    out_temp.close()
    
def parse_select_taxa():
    select_taxa = set()
    select_file = args.select_taxa_list
    
    for line in select_file:
        select_taxa.add(line.strip())
    
    select_file.close()
    return(select_taxa)    
    
if args.generate_heatmap:
    if args.significant_objects_file:
        with_sig_object()
    else:
        without_sig_object()
        
if args.find_enrichment:
    select_taxa_set = parse_select_taxa()
    taxonomy_name = args.taxonomy
    path = args.path_to_ko_files
    contrib_name = args.contrib_file
    pct_threshold = float(args.pct_threshold)
    pval_threshold = float(args.pval_threshold)

    metabolism_dict, ordered_list, ko_to_object_dict = map_kos(path)
    taxa_dict = parse_taxonomy(taxonomy_name)
    ko_dict, taxonomic_level_total, ko_tfa_dict = parse_contrib(contrib_name)
       
    second_round = apply_first_round(pct_threshold, pval_threshold)
    apply_second_round(3, pct_threshold, pval_threshold)

    for select_taxa in select_taxa_set:
        apply_second_round(3, pct_threshold, pval_threshold, select_taxa)
