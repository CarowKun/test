# coding=utf-8
import os
import re
import numpy as np
import pandas as pd
import warnings

# 不同文件的相同基因都过滤了？
os.chdir('C:/Users/25383/Desktop/Bat/gff/result_modify')
warnings.filterwarnings("ignore")


def CDS_predict(*base_seq, direction=None):
    index = []
    for z, each in enumerate(base_seq):
        index.extend(['exon' + str(z) + '|' + str(i) for i in range(len(each))])
    seq = ''.join(base_seq)
    seq = seq.upper()
    rev_seq = seq_reverse(seq)
    forword = find_length(seq)
    reverse = find_length(rev_seq)
    if len(forword) > 0 or len(reverse) > 0:
        if direction == None:
            if len(forword) >= len(reverse):
                location = re.search(forword, seq).span()
                return {'sequence': forword, 'direction': '+', 'begin': index[location[0]],
                        'end': index[location[1] - 1]}
            else:
                location = re.search(reverse, rev_seq).span()
                index = index[::-1]
                return {'sequence': reverse, 'direction': '-', 'begin': index[location[0]],
                        'end': index[location[1] - 1]}
        elif direction == '+':
            if len(forword) > 0:
                location = re.search(forword, seq).span()
                return {'sequence': forword, 'direction': '+', 'begin': index[location[0]],
                        'end': index[location[1] - 1]}
            else:
                return None
        elif direction == '-':
            if len(reverse) > 0:
                location = re.search(reverse, rev_seq).span()
                index = index[::-1]
                return {'sequence': reverse, 'direction': '-', 'begin': index[location[0]],
                        'end': index[location[1] - 1]}
            else:
                return None
    else:
        return None


def seq_reverse(base_seq):
    change_dict = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A', 'N': 'N'}
    base_seq = base_seq[::-1]
    rev_seq = [change_dict[i.upper()] for i in base_seq]
    rev_seq = ''.join(rev_seq)
    return rev_seq


gtf_wu = pd.read_csv('bat_wu.gff3', sep='\t', header=None, skiprows=[0])
gtf_yyt = pd.read_csv('bat_yyt.gff3', sep='\t', header=None, skiprows=[0])
gtf_zyx = pd.read_csv('bat_zyx.gff3', sep='\t', header=None, skiprows=[0])

gtf_wu = gtf_wu[~gtf_wu[0].isin(['HiC_scaffold_' + str(i) for i in range(7, 10)])]
gtf_yyt = gtf_yyt[gtf_yyt[0].isin(['HiC_scaffold_' + str(i) for i in '7 6 53 52 9 8'.split(' ')])]
gtf_zyx = gtf_zyx[gtf_zyx[0].isin(['HiC_scaffold_' + str(i) for i in '19 2 20 21 22 23 24'.split(' ')])]


# 删除GSAman的行
def mtx_gsaman_drop(mtx):
    mtx['gsaman'] = [False if 'GSAman00' in i else True for i in mtx[8]]  # 第八列存在GSAman标记为Fasle
    mtx = mtx[mtx['gsaman']]  # 提取第八列为True的列
    mtx = mtx.drop(columns=['gsaman'])  # 删除列，但原始的 mtx DataFrame 不会被修改
    return mtx
gtf_yyt = mtx_gsaman_drop(gtf_yyt)
gtf_zyx = mtx_gsaman_drop(gtf_zyx)


# 文件检验
def file_check(mtx):  # 检查第八列是否有Trans
    for each in mtx[8]:
        temp = re.search('.[Tt]rans[0-9]+', each)
        if temp == None:
            print(each)
file_check(gtf_yyt)
file_check(gtf_zyx)


# 不同文件Trans名改成不一致
def trans_num_move(trans, move):
    match = re.search('.[Tt]rans[0-9]+', trans).group()
    change = '.Trans' + str(int(match[6:]) + move)
    trans = trans.replace(match, change)
    return trans


gtf_yyt[8] = gtf_yyt[8].apply(trans_num_move, args=(100,))  # args传递给move；gtf_yyt[8]传递给trans；结果是Trans1——Trans101？？
gtf_zyx[8] = gtf_zyx[8].apply(trans_num_move, args=(200,))

gtf_mtx = pd.concat([gtf_wu.iloc[1:10,:], gtf_yyt.iloc[1:10,:], gtf_zyx.iloc[1:10,:]])
gtf_mtx[8] = gtf_mtx[8].str.replace('orf', 'ORF')
gtf_mtx.index = range(gtf_mtx.shape[0])


# 获取唯一gene id
genes = []
for i in gtf_mtx[8]:
    temp = re.search('=[A-Za-z0-9-\\\/]+.[Tt]rans[0-9]+', i).group()[1:]
    genes.append(temp)
gtf_mtx['genes'] = [i.split('.')[0] for i in genes]
check = gtf_mtx[[0, 'genes']]
check = check.drop_duplicates()
# print(len(check)) # 17625
# print(len(genes)) # 2080636
check = check['genes'].value_counts()
check = check[check > 1]
if check.shape[0] > 0:
    print(check)
    print('Bug')


# 更改CHUN等名
def rename_def(name1):
    for each in rename1.keys():
        name1 = name1.replace(each, rename1[each])
    return name1
record = []
for i in gtf_mtx['genes']:
    if 'ORF' in i:
        record.append(i)
record = sorted(set(record))
record3 = [i for i in record if i.startswith('CUNH')]
rename3 = dict([[i, i.replace('ORF', 'orf')] for i in record3])
record = [i for i in record if not i.startswith('CUNH')]
record1 = [i for i in record if 'H' in i]
rename1 = {}
for i in record1:
    rename = re.search('C[0-9]+', i).group()
    rename = i.replace(rename, 'CUN')
    rename = rename.replace('ORF', 'orf')
    rename1[i] = rename
record2 = [i for i in record if not 'H' in i]
# record2.remove('MORF4L2')
# record2.remove('MORF4L1')
rename2 = {}
for i in record2:
    rename = 'CUNH' + i[1:]
    rename = rename.replace('ORF', 'orf')
    rename2[i] = rename
rename1.update(rename2)
rename1.update(rename3)
gtf_mtx[8] = gtf_mtx[8].map(rename_def)  # map(rename_def)将rename_def函数映射到该列的每个元素上


# 更新过ORF所以"gene"列也要更新一遍
genes = []
for i in gtf_mtx[8]:
    temp = re.search('=[A-Za-z0-9-\\\/]+.[Tt]rans[0-9]+', i).group()[1:]
    genes.append(temp)
gtf_mtx['genes'] = [i.split('.')[0] for i in genes]


# 将第八列改为'gene_id "{gene}"; transcript_id "{trans}"的格式
genes = []
transes = []
for index in gtf_mtx.index:
    rename = gtf_mtx.loc[index, 8]
    trans = re.search('=[A-Za-z0-9-\\\/]+.[Tt]rans[0-9]+', rename).group()[1:]
    gene = trans.split('.')[0]
    genes.append(gene)
    transes.append(trans)
    new = f'gene_id "{gene}"; transcript_id "{trans}";'
    gtf_mtx.loc[index, 8] = new


# 添加额外信息
gtf_mtx['genes'] = genes
gtf_mtx['trans'] = transes
gtf_mtx[3] = gtf_mtx[3].astype(int)
gtf_mtx[4] = gtf_mtx[4].astype(int)
gtf_mtx['length'] = gtf_mtx[4] - gtf_mtx[3]
gtf_mtx[2] = [i if i != 'mRNA' else 'transcript' for i in gtf_mtx[2]]


# 构造scaffold:ATCG的字典
index_sort = []
fafile = open('../../fa/test.fa')
for line in fafile:
    if line.startswith('>'):
        index_sort.append(line[1:-1])
fafile.close()
# print(index_sort)  #['HiC_scaffold_1', 'HiC_scaffold_2', 'HiC_scaffold_3', 'HiC_scaffold_4']
fa = open('../../fa/test.fa').read()
fadict = {}
fa = fa.split('>')[
     1:]  # 取舍弃第一个是因为['', 'HiC_scaffold_1\nATCG\n', 'HiC_scaffold_2\nTCGA\n', 'HiC_scaffold_3\nCGAT\n', 'HiC_scaffold_4\nGATC']
# print(fa)
for z, line in enumerate(fa):
    line = line.split('\n', 1)
    fadict[line[0].split(' ')[0]] = line[1].replace('\n', '')
# print(fadict)  #{'HiC_scaffold_1': 'ATCG', 'HiC_scaffold_2': 'TCGA', 'HiC_scaffold_3': 'CGAT', 'HiC_scaffold_4': 'GATC'}


# 按照染色体号、基因位置排序
gtf_mtx[0] = pd.Categorical(gtf_mtx[0], categories=index_sort, ordered=True)
gtf_mtx = gtf_mtx.sort_values(by=[0, 3])
genes = gtf_mtx['genes'].drop_duplicates().tolist()


# 判断收尾想接的情况
wrong_check = pd.DataFrame()
exon = gtf_mtx[gtf_mtx[2] == 'exon']
for scaf in index_sort:
    s_sub = exon[exon[0] == scaf]
    # print(s_sub)
    for direct in ['+', '-']:
        sub = s_sub[s_sub[6] == direct]
        sub = sub.sort_values(by=[3, 4])
        if sub.shape[0] < 2:
            continue
        # print("sub=",sub)
        sub1 = sub[1:]
        sub0 = sub[:-1]
        # print("sub1=",sub1)
        # print("sub0=", sub0)

        check = [sub1[3].tolist()[i] - sub0[4].tolist()[i] for i in range(sub1.shape[0])]
        # print(check)
        sub0['gene_check'] = sub1['genes'].tolist()
        sub0['check'] = check
        # print(sub0)
        sub0 = sub0[sub0['check'] <= 0]
        sub0 = sub0[sub0['genes'] != sub0['gene_check']]
        # print(sub0)
        wrong = sub0[['genes', 'gene_check']]
        sub0['wrong'] = ['<->'.join(sorted(wrong.loc[i])) for i in wrong.index]
        wrong_check = pd.concat([wrong_check, sub0])
wrong_check.to_excel('wrong_check.xlsx', index=False)



mtx_new = pd.DataFrame()
for z, gene in enumerate(genes):
    # print(z, gene)
    sub = gtf_mtx[gtf_mtx['genes'] == gene]
    rename_mtx = pd.DataFrame(columns=['ex_len', 'ex_begin', 'ex_end', 'cds_len', 'cds_str', 'ex_str', 'direction'])
    Chr = sub.iloc[0, 0]

    new_sub = []
    trans_num = len(sub['trans'].unique())
    print(sub['trans'])

    for trans in sub['trans'].unique():
        cds_predict = False
        trans_sub = sub[sub['trans'] == trans]
        trans_sub = trans_sub.sort_values(by=3)
        exon = trans_sub[trans_sub[2] == 'exon']
        cds = trans_sub[trans_sub[2] == 'CDS']
        mrna = trans_sub[trans_sub[2] == 'transcript']
        print("mrna=",mrna)
        print("exon=",exon)

        assert mrna.iloc[0, 3] == exon.iloc[0, 3]
        assert mrna.iloc[0, 4] == exon.iloc[-1, 4]
        if cds.shape[0] > 1:
            assert exon.iloc[0, 3] <= cds.iloc[0, 3]
            assert exon.iloc[-1, 4] >= cds.iloc[0, 4]
            direction = cds.iloc[0, 6]
            if direction == '+':
                first = cds.iloc[0, 3]
                end = cds.iloc[-1, 4]
                first = fadict[Chr][first - 1:first + 2]
                end = fadict[Chr][end - 3:end]
            else:
                first = cds.iloc[-1, 4]
                end = cds.iloc[0, 3]
                first = seq_reverse(fadict[Chr][first - 3:first])
                end = seq_reverse(fadict[Chr][end - 1:end + 2])
            if first != 'ATG' or end not in ['TAA', 'TGA', 'TAG']:
                cds_predict = True
                if trans_num == 1:
                    direction = None
        else:
            cds_predict = True
            direction = None
        if cds_predict:
            exon_sequence = [fadict[Chr][exon.iloc[i, 3] - 1:exon.iloc[i, 4]] for i in range(exon.shape[0])]
            new_cds = CDS_predict(*exon_sequence, direction=direction)
            if new_cds == None:
                continue
            if new_cds['direction'] == '+':
                cds = exon[int(new_cds['begin'].split('|')[0][4:]):int(new_cds['end'].split('|')[0][4:]) + 1].copy()
                cds_copy = cds.copy()
                cds.iloc[0, 3] = cds_copy.iloc[0, 3] + int(new_cds['begin'].split('|')[1])
                cds.iloc[-1, 4] = cds_copy.iloc[-1, 3] + int(new_cds['end'].split('|')[1])
                remainder = 0
                for i in range(cds.shape[0]):
                    cds.iloc[i, 7] = str((3 - remainder) % 3)
                    remainder = (remainder + (cds.iloc[i, 4] - cds.iloc[i, 3] + 1) % 3) % 3
            else:
                cds = exon[int(new_cds['end'].split('|')[0][4:]):int(new_cds['begin'].split('|')[0][4:]) + 1].copy()
                cds_copy = cds.copy()
                cds.iloc[0, 3] = cds_copy.iloc[0, 3] + int(new_cds['end'].split('|')[1])
                cds.iloc[-1, 4] = cds_copy.iloc[-1, 3] + int(new_cds['begin'].split('|')[1])
                remainder = 0
                for i in range(-1, -1 * (cds.shape[0] + 1), -1):
                    cds.iloc[i, 7] = str((3 - remainder) % 3)
                    remainder = (remainder + (cds.iloc[i, 4] - cds.iloc[i, 3] + 1) % 3) % 3
            cds[2] = 'CDS'
            cds['length'] = cds[4] - cds[3]
            trans_sub = trans_sub[trans_sub[2] != 'CDS']
            trans_sub = pd.concat([trans_sub, cds])
            if direction == None:
                direction = new_cds['direction']
        trans_sub[6] = direction
        cds['str'] = cds[3].astype(str) + '-' + cds[4].astype(str)
        exon['str'] = exon[3].astype(str) + '-' + exon[4].astype(str)
        new_sub.append(trans_sub)
        rename_mtx.loc[trans, 'direction'] = direction
        rename_mtx.loc[trans, 'ex_len'] = np.sum(exon['length'])
        rename_mtx.loc[trans, 'ex_begin'] = exon.iloc[0, 3]
        rename_mtx.loc[trans, 'ex_end'] = exon.iloc[-1, 4]
        rename_mtx.loc[trans, 'cds_str'] = '|'.join(cds['str'])
        rename_mtx.loc[trans, 'cds_len'] = np.sum(cds['length'])
        ex_str = '|'.join(exon['str'])
        if ex_str.count('|') == 0:
            ex_str = '-'
        else:
            ex_str = ex_str.split('-', 1)[1].rsplit('-', 1)[0]
        rename_mtx.loc[trans, 'ex_str'] = ex_str
    sub = pd.concat(new_sub)
    direction = rename_mtx.loc[rename_mtx['cds_len'] == np.max(rename_mtx['cds_len']), 'direction']
    direction = list(set(direction))
    assert len(direction) == 1
    direction = direction[0]
    rename_mtx = rename_mtx[rename_mtx['direction'] == direction]
    rename_mtx = rename_mtx.sort_values(by='ex_len')
    rename_mtx = rename_mtx.drop_duplicates(subset=['cds_str', 'ex_str'])
    sub = sub[sub['trans'].isin(rename_mtx.index)]
    gene_ser = sub.iloc[[0], :].copy()
    gene_ser.iloc[0, 2] = 'gene'
    gene_ser.iloc[0, 3] = np.min(sub[3])
    gene_ser.iloc[0, 4] = np.max(sub[4])
    gene_ser.iloc[0, 8] = f'gene_id "{gene}";'
    mtx_new = pd.concat([mtx_new, gene_ser])
    rename_mtx = rename_mtx.sort_values(by=['ex_begin', 'ex_end'])
    rename_mtx['new'] = [gene + '.Trans' + str(i) for i in range(1, rename_mtx.shape[0] + 1)]
    sub.index = range(sub.shape[0])
    for each in sub.index:
        sub.loc[each, 8] = sub.loc[each, 8].replace(sub.loc[each, 'trans'], rename_mtx['new'][sub.loc[each, 'trans']])
        sub.loc[each, 'sort'] = int(rename_mtx['new'][sub.loc[each, 'trans']].split('.')[1][5:])
    sub[2] = pd.Categorical(sub[2], categories=['transcript', 'exon', 'CDS'], ordered=True)
    sub = sub.sort_values(by=['sort', 2])
    sub[2] = sub[2].tolist()
    del sub['sort']
    mtx_new = pd.concat([mtx_new, sub])


