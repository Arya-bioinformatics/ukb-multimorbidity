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

print('------------------ compare ukb and blair -----------------')
print('ukb and cell common complex diseases: ' + str(len(ukb_complex_common_disease)))
print('ukb and cell common mendelian diseases: ' + str(len(ukb_mendelian_common_disease)))
print('involved ukb multimorbidity: ' + str(len(ukb_multimorbidity)))
print('involved cell multimorbidity: ' + str(len(cell_multimorbidity)))
print('intersection: ' + str(len(ukb_multimorbidity & cell_multimorbidity)))


a = len(ukb_multimorbidity & cell_multimorbidity)
b = len(cell_multimorbidity) - a
c = len(ukb_multimorbidity) - a
d = len(ukb_complex_common_disease) * len(ukb_mendelian_common_disease) - a - b - c
[odd, p] = fisher_exact([[a,b], [c,d]])
print('-- overlap significance --')
print([odd, p])



print('ok')