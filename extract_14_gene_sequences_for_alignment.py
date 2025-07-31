# 程序功能：根据菌株ID按顺序串联14个基因的序列

gene_path_list_file = input("enter the 14 genes' path list file: ")
strain_ID_list = input('enter the strain ID list file: ')
output = input('enter the integrate_genes file: ')

list_14_genes = [
    'atp6', 'atp8', 'atp9', 'cob', 'cox1', 'cox2', 'cox3',
    'nad1', 'nad2', 'nad3', 'nad4', 'nad4L', 'nad5', 'nad6',
]

with open(strain_ID_list, 'r') as OA, open(output, 'w') as OD:
    OA_line = OA.readlines()
    for i in range(0, len(OA_line), 1):
        ID_name = OA_line[i].split('\n')[0]
        OD.write('>' + ID_name + '\n')
        for j in range(0, len(list_14_genes), 1):
            gene_ID = list_14_genes[j]
            with open(gene_path_list_file, 'r') as OB:
                OB.seek(0)
                OB_line = OB.readlines()
                for k in range(0, len(OB_line), 1):
                    path = OB_line[k].split('\n')[0]
                    with open(path, 'r') as OC:
                        OC_line = OC.readlines()
                        for m in range(0, len(OC_line), 1):
                            if f'>{ID_name}|{gene_ID}' in OC_line[m]:
                                for n in range(m+1, len(OC_line), 1):
                                    if '>' not in OC_line[n]:
                                        sequence = OC_line[n].split('\n')[0]
                                        OD.write(sequence)
                                    else:
                                        break
        OD.write('\n')
