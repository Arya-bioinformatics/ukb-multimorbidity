import pandas as pd


for flag in ['with', 'without']:
    print('\n' + flag + '\n')
    if flag == 'with':
        path1 = '../overlap/multimorbidity_snp.txt'
    if flag == 'without':
        path1 = '../overlap/multimorbidity_snp_rmhla.txt'

    multimorbidity_snp = dict()
    snp_multimorbidity_count = dict()
    disease_code_name = dict()
    with open(path1, 'r') as infile:
        for i, line in enumerate(infile):
            if i < 1:
                continue
            str1 = line.strip('\r\n')
            str2 = str1.split('\t')
            set1 = set(str2[4].split(';'))
            disease_code_name[str2[0]] = str2[2]
            disease_code_name[str2[1]] = str2[3]
            multimorbidity_snp[(str2[0], str2[1])] = set1
            for each in set1:
                if each not in snp_multimorbidity_count:
                    snp_multimorbidity_count[each] = 0
                snp_multimorbidity_count[each] += 1
        infile.close()

    list1 = list()
    for each in snp_multimorbidity_count:
        list1.append([each, snp_multimorbidity_count[each]])
    df = pd.DataFrame(list1, columns=['snp', 'count'])
    df = df.sort_values(by='count', ascending=False)

    top_snp = set()
    list1 = df.values.tolist()
    for each in list1:
        top_snp.add(each[0])
        if len(top_snp) == 10:
            break

    involved_multimorbidity = set()
    for each in multimorbidity_snp:
        name1 = disease_code_name[each[0]]
        name2 = disease_code_name[each[1]]
        if len(top_snp & multimorbidity_snp[each]) != 0:
            involved_multimorbidity.add(name1 + ' ------- ' + name2)

    for each in involved_multimorbidity:
        print(each)
    print('total snp interpreted multimorbidity: ' + str(len(multimorbidity_snp)))
    print('total top snp interpreted multimorbidity: ' + str(len(involved_multimorbidity)))
    print('ratio: ' + str(len(involved_multimorbidity)/len(multimorbidity_snp)))