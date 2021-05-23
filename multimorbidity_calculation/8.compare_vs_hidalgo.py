import pandas as pd
from itertools import product, combinations
from scipy.stats import fisher_exact
from statsmodels.stats import multitest



disease_merged = dict()
disease_description = dict()
df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    set1 = set(each[0].split(';'))
    disease_description[each[0]] = each[1]
    for each1 in set1:
        disease_merged[each1] = each[0]

df = pd.read_csv('../phenotype/multimorbidity_filter.csv')
ukb_disease = set(df['disease1']) | set(df['disease2'])


df1 = pd.read_table('../phenotype/Hidalgo/AllNet3.tsv',
                   dtype={'#icd9-1':str, 'icd9-2':str})
hidalgo_disease = set(df1['#icd9-1']) | set(df1['icd9-2'])
print('all disease: ' + str(len(hidalgo_disease)))
print('all disease pairs: ' + str(df1.shape[0]))




umls_icd10_map = dict()
with open('../name_map/2018AB-full/ICD10/2018AB/META/MRCONSO.RRF', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('|')
        umls_id = str2[0]
        icd10_id = str2[13]
        if umls_id not in umls_icd10_map:
            umls_icd10_map[umls_id] = set()
        umls_icd10_map[umls_id].add(icd10_id)
    infile.close()


icd9_icd10_map = dict()
with open('../name_map/2018AB-full/ICD9CM/2018AB/META/MRCONSO.RRF', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('|')
        umls_id = str2[0]
        icd9_id = str2[10]
        icd9_id = ''.join(['0'] * (3-len(icd9_id.split('.')[0]))) + icd9_id
        try:
            float(icd9_id)
        except:
            continue
        if umls_id in umls_icd10_map:
            icd10_set = umls_icd10_map[umls_id]
            for each in icd10_set:
                if '-' in each:
                    continue
                if icd9_id not in icd9_icd10_map:
                    icd9_icd10_map[icd9_id] = set()
                icd9_icd10_map[icd9_id].add(each)
    infile.close()

# dict1 is 3-digit icd9 map to level 2 icd10
# dict2 is 3-digit icd9 map to under level 2 icd10
# dict3 is 3-digit icd9 with their children map to level 2 icd10
# dict4 is 3-digit icd9 with their children map to under level 2 icd10
dict1 = dict()
dict2 = dict()
dict3 = dict()
dict4 = dict()
for each in hidalgo_disease:
    if each in icd9_icd10_map:
        set1 = icd9_icd10_map[each]
        set2 = set()
        set3 = set()
        for each1 in set1:
            if '-' in each1:
                continue
            if len(each1) == 3:
                set2.add(each1)
            elif len(each1) > 3:
                set3.add(each1[:3])
        if len(set2) != 0:
            dict1[each] = set2
        if len(set3) != 0:
            dict2[each] = set3
    set4 = set()
    set5 = set()
    for each1 in icd9_icd10_map:
        if each + '.' in each1:
            set1 = icd9_icd10_map[each1]
            for each2 in set1:
                if '-' in each2:
                    continue
                if len(each2) == 3:
                    set4.add(each2)
                elif len(each2) > 3:
                    set5.add(each2[:3])
    if len(set4) != 0:
        dict3[each] = set4
    if len(set5) != 0:
        dict4[each] = set5


print('---------------------------------')
print('hidalgo 3-digit icd9 map to level 2 icd10:', len(dict1), '-', len(set.union(*dict1.values())))
print('hidalgo 3-digit icd9 map to under level 2 icd10:', len(dict2), '-', len(set.union(*dict2.values())))
print('hidalgo 3-digit icd9 with their children map to level 2 icd10:', len(dict3), '-', len(set.union(*dict3.values())))
print('hidalgo 3-digit icd9 with their children map to under level 2 icd10:', len(dict4), '-', len(set.union(*dict4.values())))

print('---------------------------------')
hidalgo_disease_mapped = set(dict1.keys()) | set(dict2.keys()) | set(dict3.keys()) | set(dict4.keys())
print('hidalgo disease can map to icd10:', len(hidalgo_disease_mapped))
print('hidalgo disease with no icd10 mapped information:', len(hidalgo_disease) - len(hidalgo_disease_mapped))


icd9_icd10_childMap = dict()
for each in hidalgo_disease_mapped:
    icd9_icd10_childMap[each] = set()
    if each in dict1:
        icd9_icd10_childMap[each] |= dict1[each]
    if each in dict2:
        icd9_icd10_childMap[each] |= dict2[each]
    if each in dict3:
        icd9_icd10_childMap[each] |= dict3[each]
    if each in dict4:
        icd9_icd10_childMap[each] |= dict4[each]


print('---------------------------------')
print('hidalgo icd9-icd10 child merge:', len(icd9_icd10_childMap), '-', len(set.union(*icd9_icd10_childMap.values())))


hidalgo_disease_inUKB = dict()
ukb_disease_involved = set()
for each in icd9_icd10_childMap:
    set1 = icd9_icd10_childMap[each]
    set2 = set()
    for each1 in set1:
        if each1 in disease_merged:
            set2.add(disease_merged[each1])
            ukb_disease_involved.add(disease_merged[each1])
    if len(set2) != 0:
        hidalgo_disease_inUKB[each] = set2


ukb_disease_involved = set.union(*hidalgo_disease_inUKB.values())

# ------------ compare hidalgo and ukb (primary p values) -------------- #
hidalgo_multimorbidity = set()
df_pval = df1[(df1['RR']>1) & (df1['pval']<0.05)]
list1 = df_pval.values.tolist()
for each in list1:
    icd9_1, icd9_2 = each[0:2]
    if (icd9_1 in hidalgo_disease_inUKB) & (icd9_2 in hidalgo_disease_inUKB):
        set1 = hidalgo_disease_inUKB[icd9_1]
        set2 = hidalgo_disease_inUKB[icd9_2]
        for icd10_1, icd10_2 in product(set1, set2):
            if icd10_1 < icd10_2:
                hidalgo_multimorbidity.add(icd10_1 + '-' + icd10_2)
            else:
                hidalgo_multimorbidity.add(icd10_2 + '-' + icd10_1)

ukb_multimorbidity = set()
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[:2]
    if (code1 in ukb_disease_involved) & (code2 in ukb_disease_involved):
        ukb_multimorbidity.add(code1 + '-' + code2)

print('------------- compare hidalgo and ukb (primary p values) ---------------')
print('hidalgo disease map to ukb disease:', len(hidalgo_disease_inUKB), '-', len(set.union(*hidalgo_disease_inUKB.values())))
print('involved hidalgo multimorbidity: ' + str(len(hidalgo_multimorbidity)))
print('involved ukb multimorbidity: ' + str(len(ukb_multimorbidity)))
print('hidalgo and ukb multimorbidity overlap: ' + str(len(hidalgo_multimorbidity & ukb_multimorbidity)))


a = len(ukb_multimorbidity & hidalgo_multimorbidity)
b = len(hidalgo_multimorbidity) - a
c = len(ukb_multimorbidity) - a
d = len(set(combinations(ukb_disease_involved, 2))) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])



# ------------ compare hidalgo and ukb (fdr correction) -------------- #
hidalgo_multimorbidity = set()
df2 = df1[df1['RR']>1]
df2['q'] = multitest.fdrcorrection(df2['pval'])[1]
df_fdr = df2[df2['q']<0.05]


list1 = df_fdr.values.tolist()
for each in list1:
    icd9_1, icd9_2 = each[0:2]
    if (icd9_1 in hidalgo_disease_inUKB) & (icd9_2 in hidalgo_disease_inUKB):
        set1 = hidalgo_disease_inUKB[icd9_1]
        set2 = hidalgo_disease_inUKB[icd9_2]
        for icd10_1, icd10_2 in product(set1, set2):
            if icd10_1 < icd10_2:
                hidalgo_multimorbidity.add(icd10_1 + '-' + icd10_2)
            else:
                hidalgo_multimorbidity.add(icd10_2 + '-' + icd10_1)


ukb_multimorbidity = set()
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[:2]
    if (code1 in ukb_disease_involved) & (code2 in ukb_disease_involved):
        ukb_multimorbidity.add(code1 + '-' + code2)

print('------------ compare hidalgo and ukb (fdr correction) --------------')
print('hidalgo disease map to ukb disease:', len(hidalgo_disease_inUKB), '-', len(set.union(*hidalgo_disease_inUKB.values())))
print('involved hidalgo multimorbidity: ' + str(len(hidalgo_multimorbidity)))
print('involved ukb multimorbidity: ' + str(len(ukb_multimorbidity)))
print('hidalgo and ukb multimorbidity overlap: ' + str(len(hidalgo_multimorbidity & ukb_multimorbidity)))


a = len(ukb_multimorbidity & hidalgo_multimorbidity)
b = len(hidalgo_multimorbidity) - a
c = len(ukb_multimorbidity) - a
d = len(set(combinations(ukb_disease_involved, 2))) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])








print('ok')