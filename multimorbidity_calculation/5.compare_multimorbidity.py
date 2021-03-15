import pandas as pd
from itertools import combinations, product
from scipy.stats import fisher_exact
from statsmodels.stats import multitest




disease_merged = dict()
disease_description = dict()
df = pd.read_excel('/../phenotype/disease_class_manual.xlsx')
list1 = df.values.tolist()
for each in list1:
    set1 = set(each[0].split(';'))
    disease_description[each[0]] = each[1]
    for each1 in set1:
        disease_merged[each1] = each[0]


df = pd.read_csv('/../phenotype/multimorbidity_filter.csv')
ukb_disease = set(df['disease1']) | set(df['disease2'])


df1 = pd.read_table('/../phenotype/multimorbidity_Hidalgo/AllNet3.tsv',
                   dtype={'#icd9-1':str, 'icd9-2':str})
hidalgo_disease = set(df1['#icd9-1']) | set(df1['icd9-2'])
print('all disease: ' + str(len(hidalgo_disease)))
print('all disease pairs: ' + str(df1.shape[0]))


icd9_icd10_map = dict()
icd9_code = set()
icd10_code = set()
with open('/../icd10_icd9cm_map.txt', 'r') as infile:
    for line in infile:
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        icd10 = str2[0]
        icd10_code.add(icd10)
        set1 = set(str2[1].split(';'))
        icd9_code |= set1
        for each in set1:
            if each not in icd9_icd10_map:
                icd9_icd10_map[each] = set()
            icd9_icd10_map[each].add(icd10)
    infile.close()

print('icd9 codes mapped to icd10 in umls: ' + str(len(icd9_code)))
print('icd10 codes mapped to icd9 in umls: ' + str(len(icd10_code)))

icd9_icd10_3digital = dict()
for each in icd9_icd10_map:
    if len(each) != 3:
        continue
    set1 = icd9_icd10_map[each]
    set2 = set()
    for each1 in set1:
        if len(each1) == 3:
            set2.add(each1)
    if len(set2) != 0:
        icd9_icd10_3digital[each] = ''.join(set2)
        if len(set2) > 1:
            print('number of icd10 codes mapped by icd9: ' + str(len(set2)))

hidalgo_icd10_disease = set()
set1 = set()
set2 = set()
for each in hidalgo_disease:
    if each in icd9_icd10_3digital:
        hidalgo_icd10_disease.add(icd9_icd10_3digital[each])
    elif each in icd9_icd10_map:
        set1.add(each)
    else:
        set2.add(each)

print('hidalgo diseases mapped to level 2 icd10: ' + str(len(hidalgo_icd10_disease)))
print('hidalgo diseases mapped to other level icd10: ' + str(len(set1)))
print('hidalgo diseases can not be mapped to icd10: ' + str(len(set2)))

with open('a.csv', 'w+') as outfile:
    for each in set2:
        outfile.write(each + '\n')
    outfile.close()


hidalgo_merged_disease = set()
for each in hidalgo_icd10_disease:
    if each in disease_merged:
        hidalgo_merged_disease.add(disease_merged[each])
    else:
        hidalgo_merged_disease.add(each)


# ------------ compare hidalgo and ukb (primary p values) -------------- #
hidalgo_ukb_common_disease = hidalgo_merged_disease & ukb_disease

hidalgo_multimorbidity = set()
df_pval = df1[(df1['RR']>1) & (df1['pval']<0.05)]
list1 = df_pval.values.tolist()
for each in list1:
    icd9_1, icd9_2 = each[0:2]
    if (icd9_1 in icd9_icd10_3digital) & (icd9_2 in icd9_icd10_3digital):
        icd10_1 = icd9_icd10_3digital[icd9_1]
        icd10_2 = icd9_icd10_3digital[icd9_2]
        if icd10_1 in disease_merged:
            icd10_1 = disease_merged[icd10_1]
        if icd10_2 in disease_merged:
            icd10_2 = disease_merged[icd10_2]
        if (icd10_1 in hidalgo_ukb_common_disease) & (icd10_2 in hidalgo_ukb_common_disease):
            if icd10_1 < icd10_2:
                hidalgo_multimorbidity.add(icd10_1 + '-' + icd10_2)
            else:
                hidalgo_multimorbidity.add(icd10_2 + '-' + icd10_1)


ukb_multimorbidity = set()
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[:2]
    if (code1 in hidalgo_ukb_common_disease) & (code2 in hidalgo_ukb_common_disease):
        ukb_multimorbidity.add(code1 + '-' + code2)

print('------------- compare hidalgo and ukb (primary p values) ---------------')
print('ukb and hidalgo common diseases: ' + str(len(hidalgo_ukb_common_disease)))
print('involved hidalgo multimorbidity: ' + str(len(hidalgo_multimorbidity)))
print('involved ukb multimorbidity: ' + str(len(ukb_multimorbidity)))
print('hidalgo and ukb multimorbidity overlap: ' + str(len(hidalgo_multimorbidity & ukb_multimorbidity)))


a = len(ukb_multimorbidity & hidalgo_multimorbidity)
b = len(hidalgo_multimorbidity) - a
c = len(ukb_multimorbidity) - a
d = len(set(combinations(hidalgo_ukb_common_disease, 2))) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])




# ------------ compare hidalgo and ukb (fdr correction) -------------- #
hidalgo_ukb_common_disease = hidalgo_merged_disease & ukb_disease

hidalgo_multimorbidity = set()
df2 = df1[df1['RR']>1]
df2['q'] = multitest.fdrcorrection(df2['pval'])[1]
df_fdr = df2[df2['q']<0.05]
list1 = df_fdr.values.tolist()
for each in list1:
    icd9_1, icd9_2 = each[0:2]
    if (icd9_1 in icd9_icd10_3digital) & (icd9_2 in icd9_icd10_3digital):
        icd10_1 = icd9_icd10_3digital[icd9_1]
        icd10_2 = icd9_icd10_3digital[icd9_2]
        if icd10_1 in disease_merged:
            icd10_1 = disease_merged[icd10_1]
        if icd10_2 in disease_merged:
            icd10_2 = disease_merged[icd10_2]
        if (icd10_1 in hidalgo_ukb_common_disease) & (icd10_2 in hidalgo_ukb_common_disease):
            if icd10_1 < icd10_2:
                hidalgo_multimorbidity.add(icd10_1 + '-' + icd10_2)
            else:
                hidalgo_multimorbidity.add(icd10_2 + '-' + icd10_1)


ukb_multimorbidity = set()
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[:2]
    if (code1 in hidalgo_ukb_common_disease) & (code2 in hidalgo_ukb_common_disease):
        ukb_multimorbidity.add(code1 + '-' + code2)

print('------------ compare hidalgo and ukb (fdr correction) --------------')
print('ukb and hidalgo common diseases: ' + str(len(hidalgo_ukb_common_disease)))
print('involved hidalgo multimorbidity: ' + str(len(hidalgo_multimorbidity)))
print('involved ukb multimorbidity: ' + str(len(ukb_multimorbidity)))
print('hidalgo and ukb multimorbidity overlap: ' + str(len(hidalgo_multimorbidity & ukb_multimorbidity)))


a = len(ukb_multimorbidity & hidalgo_multimorbidity)
b = len(hidalgo_multimorbidity) - a
c = len(ukb_multimorbidity) - a
d = len(set(combinations(hidalgo_ukb_common_disease, 2))) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])



# ----------------- compare hidalgo and ukb (bonferroni correction) ---------------- #
hidalgo_multimorbidity_bf = set()
df2 = df1[df1['RR']>1]
df_bf = df2[df2['pval']<(0.05/df2.shape[0])]
list1 = df_bf.values.tolist()
for each in list1:
    icd9_1, icd9_2 = each[0:2]
    if (icd9_1 in icd9_icd10_3digital) & (icd9_2 in icd9_icd10_3digital):
        icd10_1 = icd9_icd10_3digital[icd9_1]
        icd10_2 = icd9_icd10_3digital[icd9_2]
        if icd10_1 in disease_merged:
            icd10_1 = disease_merged[icd10_1]
        if icd10_2 in disease_merged:
            icd10_2 = disease_merged[icd10_2]
        if (icd10_1 in hidalgo_ukb_common_disease) & (icd10_2 in hidalgo_ukb_common_disease):
            if icd10_1 < icd10_2:
                hidalgo_multimorbidity_bf.add(icd10_1 + '-' + icd10_2)
            else:
                hidalgo_multimorbidity_bf.add(icd10_2 + '-' + icd10_1)

print('\nhidalgo and ukb multimorbidity overlap (bonferroni correction): '
      + str(len(hidalgo_multimorbidity_bf & ukb_multimorbidity)) + '\n')



# ------------------------ compare jensen and ukb ----------------------------- #
jensen_disease = set()
df3 = pd.read_excel('/../phenotype/multimorbidity_Jensen/ncomms5022-s2.xlsx')
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    jensen_disease.add(code1)
    jensen_disease.add(code2)

ukb_jensen_common_disease = ukb_disease & jensen_disease

ukb_multimorbidity = set()
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[0], each[1]
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    ukb_multimorbidity.add(code1 + '-' + code2)

jensen_multimorbidity = set()
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    if code1 < code2:
        jensen_multimorbidity.add(code1 + '-' + code2)
    else:
        jensen_multimorbidity.add(code2 + '-' + code1)

print('------------------ compare ukb and jensen -----------------')
print('ukb and jensen common diseases: ' + str(len(ukb_jensen_common_disease)))
print('involved ukb multimorbidity: ' + str(len(ukb_multimorbidity)))
print('involved jensen multimorbidity: ' + str(len(jensen_multimorbidity)))
print('intersaction: ' + str(len(ukb_multimorbidity & jensen_multimorbidity)))

a = len(ukb_multimorbidity & jensen_multimorbidity)
b = len(jensen_multimorbidity) - a
c = len(ukb_multimorbidity) - a
d = len(set(combinations(ukb_jensen_common_disease, 2))) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])




# ----------------- compare directional ukb and jensen ------------------ #
ukb_jensen_common_disease = ukb_disease & jensen_disease
jensen_disease = set()
df3 = pd.read_excel('/../phenotype/multimorbidity_Jensen/ncomms5022-s2.xlsx')
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    jensen_disease.add(code1)
    jensen_disease.add(code2)

df4 = pd.read_csv('/../phenotype/multimorbidity_filter_direction.csv')
ukb_directional_disease = set(df4.loc[df4['direction']!=0, 'disease1']) | set(df4.loc[df4['direction']!=0, 'disease2'])


ukb_multimorbidity = set()
list1 = df4.values.tolist()
for each in list1:
    code1, code2, direction = each[0], each[1], each[12]
    if direction == 0:
        continue
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    ukb_multimorbidity.add(code1 + '-' + code2)

jensen_multimorbidity = set()
list1 = df3.values.tolist()
for each in list1[4:]:
    code1, code2 = each[0], each[2]
    if code1 in disease_merged:
        code1 = disease_merged[code1]
    if code2 in disease_merged:
        code2 = disease_merged[code2]
    if code1 not in ukb_jensen_common_disease:
        continue
    if code2 not in ukb_jensen_common_disease:
        continue
    if code1 < code2:
        jensen_multimorbidity.add(code1 + '-' + code2)
    else:
        jensen_multimorbidity.add(code2 + '-' + code1)

print('------------------ compare directional ukb and jensen -----------------')
print('ukb and jensen common diseases: ' + str(len(ukb_jensen_common_disease)))
print('involved ukb multimorbidity: ' + str(len(ukb_multimorbidity)))
print('involved jensen multimorbidity: ' + str(len(jensen_multimorbidity)))
print('intersaction: ' + str(len(ukb_multimorbidity & jensen_multimorbidity)))

a = len(ukb_multimorbidity & jensen_multimorbidity)
b = len(jensen_multimorbidity) - a
c = len(ukb_multimorbidity) - a
d = len(set(combinations(ukb_jensen_common_disease, 2))) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])




# ----------------- compare ukb and blair ------------------ #
disease_name_code_map = dict()
complex_disease = set()
df4 = pd.read_excel('/../1-s2.0-S0092867413010246-mmc2.xls',encoding='ISO-8859-1')
list1 = df4.values.tolist()
for each in list1:
    name = str.lower(each[0])
    set1 = set(each[2].replace(' ', '').split(','))
    for each1 in set1:
        if each1[:3] in disease_merged:
            complex_disease.add(disease_merged[each1[:3]])
            if name not in disease_name_code_map:
                disease_name_code_map[name] = set()
            disease_name_code_map[name].add(disease_merged[each1[:3]])
        else:
            complex_disease.add(each1[:3])
            if name not in disease_name_code_map:
                disease_name_code_map[name] = set()
            disease_name_code_map[name].add(each1[:3])
mendelian_disease = set()
df5 = pd.read_excel('/../1-s2.0-S0092867413010246-mmc3.xls',encoding='ISO-8859-1')
list1 = df5.values.tolist()
for each in list1:
    name = str.lower(each[0])
    set1 = set(each[2].replace(' ', '').split(','))
    for each1 in set1:
        if each1 in disease_merged:
            mendelian_disease.add(disease_merged[each1[:3]])
            if name not in disease_name_code_map:
                disease_name_code_map[name] = set()
            disease_name_code_map[name].add(disease_merged[each1[:3]])
        else:
            mendelian_disease.add(each1[:3])
            if name not in disease_name_code_map:
                disease_name_code_map[name] = set()
            disease_name_code_map[name].add(each1[:3])

ukb_complex_common_disease = ukb_disease & complex_disease - (complex_disease & mendelian_disease)
ukb_mendelian_common_disease = ukb_disease & mendelian_disease - (complex_disease & mendelian_disease)

ukb_multimorbidity = set()
list1 = df.values.tolist()
for each in list1:
    code1, code2 = each[0], each[1]
    if (len(set([code1, code2]) & ukb_complex_common_disease) == 1) & (len(set([code1, code2]) & ukb_mendelian_common_disease) == 1):
        ukb_multimorbidity.add(code1 + '-' + code2)

cell_multimorbidity = set()
df6 = pd.read_excel('/../1-s2.0-S0092867413010246-mmc4.xls')
list1 = df6.values.tolist()
for each in list1[1:]:
    disease1 = str.lower(each[0])
    disease2 = str.lower(each[1])
    if disease1 not in disease_name_code_map:
        continue
    if disease2 not in disease_name_code_map:
        continue
    set1 = disease_name_code_map[disease1] & (ukb_complex_common_disease | ukb_mendelian_common_disease)
    set2 = disease_name_code_map[disease2] & (ukb_complex_common_disease | ukb_mendelian_common_disease)
    if (len(set1) == 0) | (len(set2) == 0):
        continue
    set3 = product(set1, set2)
    for code1, code2 in set3:
        if code1 < code2:
            cell_multimorbidity.add(code1 + '-' + code2)
        else:
            cell_multimorbidity.add(code2 + '-' + code1)

print('------------------ compare ukb and cell -----------------')
print('ukb and cell common complex diseases: ' + str(len(ukb_complex_common_disease)))
print('ukb and cell common mendelian diseases: ' + str(len(ukb_mendelian_common_disease)))
print('involved ukb multimorbidity: ' + str(len(ukb_multimorbidity)))
print('involved cell multimorbidity: ' + str(len(cell_multimorbidity)))
print('intersaction: ' + str(len(ukb_multimorbidity & cell_multimorbidity)))


a = len(ukb_multimorbidity & cell_multimorbidity)
b = len(cell_multimorbidity) - a
c = len(ukb_multimorbidity) - a
d = len(ukb_complex_common_disease) * len(ukb_mendelian_common_disease) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])



print('ok')