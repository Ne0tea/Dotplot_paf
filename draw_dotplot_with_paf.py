'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-05-27 19:21:44
LastEditors: Ne0tea
LastEditTime: 2024-05-31 16:54:21
'''
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import logging

def draw_dotplot_with_highlight(dataframe,yregion,xregion,query,ref,store_dir=os.getcwd()):
    # 绘制点图
    plt.figure(figsize=(10, 10))
    ax = plt.gca()

    # 添加 x 轴区域的背景色
    ax.axvspan(xregion[0], xregion[1], color='#003049', alpha=0.25)

    # 添加 y 轴区域的背景色
    ax.axhspan(yregion[0], yregion[1], color='#780000', alpha=0.25)
    for _, row in dataframe.iterrows():
        if row['strand'] == '+':
            plt.plot( [row['target_start'], row['target_end']], [row['query_start'], row['query_end']],
                    markersize=2,marker='o', linestyle='-',linewidth=2,color='black')
        elif row['strand'] == '-':
            plt.plot( [row['target_end'],row['target_start']], [row['query_start'], row['query_end']],
                    markersize=2,marker='o', linestyle='-',linewidth=2,color='black')
    # plt.title('Dot Plot of PAF Alignments')
    plt.xlabel(ref)
    plt.ylabel(query)
    plt.grid(True)
    plt.savefig(os.path.join(store_dir,query+'_vs_'+ref+'.pdf'), format='pdf')
    # plt.show()

def main(paf_file,centro_file,filter_len,store_dir):
    # 定义PAF文件的列名
    if os.path.exists(store_dir):
        pass
    else:
        os.makedirs(store_dir,exist_ok=True)
    centro_dic={}
    with open(centro_file,'r') as CRf:
        for i in CRf:
            line=i.strip().split('\t')
            centro_dic[line[0]]=[int(line[1]),int(line[2])]

    column_names = [
        "query_name", "query_length", "query_start", "query_end", "strand",
        "target_name", "target_length", "target_start", "target_end",
        "residue_matches", "alignment_block_length", "mapping_quality"
    ]

    # 读取PAF文件
    df = pd.read_table(paf_file,  header=None, names=column_names,usecols=range(12))
    df = df[df['alignment_block_length']>filter_len]
    plot_pre=df['query_name'].unique()

    for i in plot_pre:
        cur_chr=i.split('_')[0]
        cur_df = df[((df['query_name'].str.contains(cur_chr,na=False)) & (df['target_name'].str.contains(cur_chr,na=False)))]
        logging.info(i+'start plotting!')
        if len(cur_df['target_name'].unique()) > 1:
            cur_target=cur_df['target_name'].value_counts().idxmax()
        else:
            cur_target=cur_df['target_name'].unique()[0]
        # print(cur_target)
        query_centro=centro_dic[i]
        target_centro=centro_dic[cur_target]
        draw_dotplot_with_highlight(cur_df,query_centro,target_centro,i,cur_target,store_dir)
        # print(cur_df)
    # 提取比对起始和终止位置

if __name__ == "__main__":
    # paf_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\Chr06_ali_BS_NIP\test.paf'
    '''
    centro_dic format
    Chr01_CW11	17450634	17569863
    Chr01_CW13	19349079	19805987
    Chr01_CW14	18383125	18884740
    Chr01_CW15	18356684	18885175
    Chr01_CW18	19813075	20337906
    Chr01_CW22	17873425	18382707
    '''
    # centro_file=r'E:\Bio_analysis\Weedyrice\pan_weedyrice\CR_size_varition\centromere_region_70m_used.bed'
    parser = argparse.ArgumentParser(description="draw dot plot based on paf")
    parser.add_argument('-i',dest="paf", type=str, help="path of paf file",required=True)
    parser.add_argument('-r',dest="region", type=str, help="region file indicate high light region")
    parser.add_argument('-op',dest="outdir", default=os.getcwd(), type=str, help="outdir of plot file")
    parser.add_argument("-flen",dest="flen",default=2000, type=int, help="filter alignlen less then")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, \
                        format='%(asctime)s - %(levelname)s - %(message)s')
    main(args.paf,args.region,args.flen,args.outdir)
