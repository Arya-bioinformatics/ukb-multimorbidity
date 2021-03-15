import pandas as pd
from collections import Counter

# ----------------------------- find top disease --------------------------- #
top = 0.25

multimorbidity_loci = set()
multimorbidity_snp = dict()
with open('../overlap/multimorbidity_snp.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        multimorbidity_loci.add((code1, code2))
        multimorbidity_snp[(code1, code2)] = set(str2[4].split(';'))
    infile.close()

multimorbidity_gene = dict()
with open('../overlap/multimorbidity_gene.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        multimorbidity_loci.add((code1, code2))
        multimorbidity_gene[(code1, code2)] = set(str2[5].split(';'))
    infile.close()

disease_degree = dict()
for each in multimorbidity_loci:
    code1, code2 = each
    if code1 not in disease_degree:
        disease_degree[code1] = 0
    disease_degree[code1] += 1
    if code2 not in disease_degree:
        disease_degree[code2] = 0
    disease_degree[code2] += 1

list1 = list()
for each in disease_degree:
    list1.append([each, disease_degree[each]])
df = pd.DataFrame(list1, columns=['disease', 'degree'])
top_025_disease_loci = set(df.loc[df['degree']>df.shape[0]*top, 'disease'])
print('total top 25% disease loci: ' + str(len(top_025_disease_loci)))

multimorbidity_network = set()
multimorbidity_ppi = dict()
with open('../overlap/multimorbidity_ppi.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        multimorbidity_network.add((code1, code2))
        multimorbidity_ppi[(code1, code2)] = set(str2[5].split(';'))
    infile.close()

multimorbidity_pathway = dict()
with open('../overlap/multimorbidity_pathway.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        code1, code2 = str2[:2]
        multimorbidity_network.add((code1, code2))
        multimorbidity_pathway[(code1, code2)] = set(str2[4].split(';'))
    infile.close()

disease_degree = dict()
for each in multimorbidity_network:
    code1, code2 = each
    if code1 not in disease_degree:
        disease_degree[code1] = 0
    disease_degree[code1] += 1
    if code2 not in disease_degree:
        disease_degree[code2] = 0
    disease_degree[code2] += 1

list1 = list()
for each in disease_degree:
    list1.append([each, disease_degree[each]])
df = pd.DataFrame(list1, columns=['disease', 'degree'])
top_025_disease_network = set(df.loc[df['degree'] > df.shape[0] * top, 'disease'])
print('total top 25% disease network: ' + str(len(top_025_disease_network)))


top_disease_multimorbidity_loci = dict()
for each in top_025_disease_loci:
    set1 = set()
    for code1, code2 in multimorbidity_loci:
        if code1 == each:
            set1.add(code2)
        elif code2 == each:
            set1.add(code1)
    top_disease_multimorbidity_loci[each] = set1

top_disease_snp = dict()
for each in top_025_disease_loci:
    list1 = list()
    for each1 in multimorbidity_snp:
        if each in each1:
            list1 += list(multimorbidity_snp[each1])
    list2 = list()
    for each1 in set(list1):
        list2.append([each1, list1.count(each1)])
    df = pd.DataFrame(list2, columns=['snp', 'count'])
    list3 = list(set(df['count']))
    list3.sort(reverse=True)
    top_5_num = list3[:5]
    list3 = df.values.tolist()
    set1 = set()
    for each1 in list3:
        if each1[1] in top_5_num:
            set1.add(each1[0])
    top_disease_snp[each] = set1


top_disease_gene = dict()
for each in top_025_disease_loci:
    list1 = list()
    for each1 in multimorbidity_gene:
        if each in each1:
            list1 += list(multimorbidity_gene[each1])
    list2 = list()
    for each1 in set(list1):
        list2.append([each1, list1.count(each1)])
    df = pd.DataFrame(list2, columns=['gene', 'count'])
    list3 = list(set(df['count']))
    list3.sort(reverse=True)
    top_5_num = list3[:5]
    list3 = df.values.tolist()
    set1 = set()
    for each1 in list3:
        if each1[1] in top_5_num:
            set1.add(each1[0])
    top_disease_gene[each] = set1

result = list()
for each in top_025_disease_loci:
    set1 = top_disease_multimorbidity_loci[each]
    temp = list(set1)
    temp.sort()
    result.append([each, 'LG-network', len(set1), ' | '.join(temp), ', '.join(top_disease_gene[each])])


top_disease_multimorbidity_network = dict()
for each in top_025_disease_network:
    set1 = set()
    for code1, code2 in multimorbidity_network:
        if code1 == each:
            set1.add(code2)
        elif code2 == each:
            set1.add(code1)
    top_disease_multimorbidity_network[each] = set1

top_disease_ppi = dict()
for each in top_025_disease_network:
    list1 = list()
    for each1 in multimorbidity_ppi:
        if each in each1:
            list1 += list(multimorbidity_ppi[each1])
    list2 = list()
    for each1 in set(list1):
        list2.append([each1, list1.count(each1)])
    df = pd.DataFrame(list2, columns=['ppi', 'count'])
    list3 = list(set(df['count']))
    list3.sort(reverse=True)
    top_5_num = list3[:5]
    list3 = df.values.tolist()
    set1 = set()
    for each1 in list3:
        if each1[1] in top_5_num:
            set1.add(each1[0])
    top_disease_ppi[each] = set1


top_disease_pathway = dict()
for each in top_025_disease_network:
    list1 = list()
    for each1 in multimorbidity_pathway:
        if each in each1:
            list1 += list(multimorbidity_pathway[each1])
    list2 = list()
    for each1 in set(list1):
        list2.append([each1, list1.count(each1)])
    df = pd.DataFrame(list2, columns=['pathway', 'count'])
    list3 = list(set(df['count']))
    list3.sort(reverse=True)
    top_5_num = list3[:5]
    list3 = df.values.tolist()
    set1 = set()
    for each1 in list3:
        if each1[1] in top_5_num:
            set1.add(each1[0])
    top_disease_pathway[each] = set1


for each in top_025_disease_network:
    set1 = top_disease_multimorbidity_network[each]
    temp = list(set1)
    temp.sort()
    result.append([each, 'NG-network', len(set1), ' | '.join(temp), ', '.join(top_disease_ppi[each])])
    result.append([each, 'NG-network', len(set1), ' | '.join(temp), ', '.join(top_disease_pathway[each])])

df = pd.DataFrame(result, columns=['Universal hub disease', 'Network',
                                  'Number of interpretable multimorbidities', 'multimorbid diseases', 'Top 5 genetic components'])
df.to_csv('a.csv', index=False)