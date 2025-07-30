cds_or_pep_file = input('enter the cds_or_pep_file : ')
output = input('enter the core_genes_cds or core_genes_pep file : ')
strain_ID = input('enter the strain ID : ')

list_14_genes = [
    'atp6',
    'atp8',
    'atp9',
    'cob',
    'cox1',
    'cox2',
    'cox3',
    'nad1',
    'nad2',
    'nad3',
    'nad4',
    'nad4L',
    'nad5',
    'nad6',
]

OB = open(cds_or_pep_file, 'r')
OC = open(output, 'w')
OC.write('>'+strain_ID+'\n')
for i in range(0, len(list_14_genes), 1):
    gene_ID = list_14_genes[i]
    OB.seek(0)                               #程序的"灵魂"
    OB_line = OB.readline()
    while OB_line:
        if '[' + 'gene' + '=' + gene_ID + ']' in OB_line:
            #OC.write('>' + strain_ID + '|' + str(gene_ID) + '\n')
            OB_line = OB.readline()
            while '>' not in OB_line or OB_line != '':
                OC.write(OB_line)
                OB_line = OB.readline()
                if '>' in OB_line or OB_line == '':
                    break
        OB_line = OB.readline()
OB.close()
OC.close()


