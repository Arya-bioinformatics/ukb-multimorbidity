import pandas as pd
from scipy.stats import pearsonr




df = pd.read_csv('/../Disease_summary_information.csv')
disease_prevalence = df.set_index(df['ICD10'])['Number of patients'].to_dict()

disease_snp = dict()
with open('/../genome/disease_snp_sig.txt', 'r') as infile:
    for i, line in enumerate(infile):
        if i < 1:
            continue
        str1 = line.strip('\r\n')
        str2 = str1.split('\t')
        disease_snp[str2[0]] = set(str2[1].split(';'))
    infile.close()


list1 = list()
list2 = list()
for each in disease_snp:
    list1.append(len(disease_snp[each]))
    list2.append(disease_prevalence[each]/410293)


print(pearsonr(list1, list2))