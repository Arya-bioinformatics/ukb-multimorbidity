
import pandas as pd


if __name__ == '__main__':

    HLA_region = 'chr6:29691116–33054976'
    start = 29691116
    end = 33054976


    disease_class = dict()
    all_class = set()
    df = pd.read_excel('../phenotype/disease_class_manual.xlsx')
    list1 = df.values.tolist()
    for each in list1:
        code = each[0]
        c = each[2]
        disease_class[code] = c
        all_class.add(c)

    multimorbidity_snp = set()
    multimorbidity_snp_dic = dict()
    with open('../overlap/multimorbidity_snp.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            code1, code2 = str2[:2]
            c1 = disease_class[code1]
            c2 = disease_class[code2]
            snp = str2[4].split(';')
            multimorbidity_snp |= set(snp)
            multimorbidity_snp_dic[(code1, code2)] = set(snp)
        infile.close()


    snp_position = dict()
    with open('../genome/varaint_info.txt', 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            variant = str2[0]
            if variant not in multimorbidity_snp:
                continue
            chr = str2[1]
            position = str2[2]
            snp_position[variant] = [chr, position]
        infile.close()

    set1 = set()
    set2 = set()
    for each in multimorbidity_snp:
        [chr, pos] = snp_position[each]
        if chr != '6':
            continue
        set1.add(each)
        if (int(pos) >= start) & (int(pos) <= end):
            set2.add(each)



    print('----------- multimorbidity mhc snp ratio -------------')
    print(len(set1))
    print(len(set2))
    print(float(len(set2))/len(set1))
    print(len(set2)/len(multimorbidity_snp))

    print('----------- multimorbidity with snp in chr6, but not hla -------------')
    for each in multimorbidity_snp_dic:
        set3 = multimorbidity_snp_dic[each]
        set4 = set1 & set3
        set5 = set2 & set3
        if (len(set4) != 0) & (len(set5) == 0):
            print(each)

    print('----------- multimorbidity with hla snp overlap -------------')
    i = 0
    for each in multimorbidity_snp_dic:
        set3 = multimorbidity_snp_dic[each]
        set5 = set2 & set3
        if len(set5) != 0:
            i += 1
    print(i)
    print(len(multimorbidity_snp_dic))


    print('----------- top 10 multimorbidities with the pargest number of hla snps -------------')
    list1 = list()
    for each in multimorbidity_snp_dic:
        count = len(multimorbidity_snp_dic[each] & set2)
        list1.append([each[0], each[1], count])
    df = pd.DataFrame(list1, columns=['code1', 'code2', 'hla snp count'])
    df = df.sort_values(by='hla snp count', ascending=False)