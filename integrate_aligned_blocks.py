import argparse

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='Concatenate blocks files by matching IDs')

    # 添加命令行参数
    parser.add_argument('--strainID', required=True, help='one file containing IDs of all strain')
    parser.add_argument('--block_path_list', required=True, help='one file containing paths of all blocks file')
    parser.add_argument('--output', required=True, default='08.integrate_c_muscle.txt', help='output file name')

    # 解析命令行参数
    args = parser.parse_args()

    # 获取参数值
    strainID_file = args.strainID
    muscle_block_list_file = args.block_path_list
    output_file = args.output

    # 创建一个字典来存储菌株ID和对应的序列
    id_sequence_dict = {}

    # 读取菌株ID列表文件
    with open(strainID_file, 'r') as F:
        ID_content = F.readlines()
        for ID_line in ID_content:
            ID = ID_line.strip()
            sequence = ''
            # 读取包含Muscle比对结果路径的文件
            with open(muscle_block_list_file, 'r') as OB:
                list_content = OB.readlines()
                for list_line in list_content:
                    # 获取每一个one_block_muscle的文件路径
                    one_block_file = list_line.strip()
                    # 读取这个路径的文件
                    with open(one_block_file, 'r') as OC:
                        seq_file = OC.readlines()
                        found = False
                        for k, line in enumerate(seq_file):
                            if '>' in line:
                                ID_2 = line.split(' ')[0].split('>')[1].strip()
                                if ID_2 == ID:
                                    found = True
                                    continue
                            if found:
                                if '>' in line:
                                    break
                                sequence += line.strip()
            id_sequence_dict[ID] = sequence

    # 将结果写入输出文件
    with open(output_file, 'w') as OD:
        for ID, seq in id_sequence_dict.items():
            OD.write(f'>{ID}\n')
            OD.write(f'{seq}\n')

if __name__ == "__main__":
    main()

