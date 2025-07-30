# 引入模块，可以在python系统下运行linux命令
import os
import subprocess
import argparse
from collections import defaultdict

# 定义函数：把一个文件的内容追加到另一个文件中
def copy_file(file_1, file_2):
    # 读取文件1
    AA = open(file_1, 'r')
    # 按行读取
    file_1_content = AA.readlines()
    # 打开文件2，追加
    BB = open(file_2, 'a+')
    # 遍历列表
    for x in range(0, len(file_1_content), 1):
        # 将文件1的每行追加到文件2中
        BB.write(file_1_content[x])
    AA.close()
    BB.close()

# 定义函数：读取FASTA格式文件获得特定ID的序列长度
# 定义函数ref_genome_size，所需参数file_path（文件路径）；reference_name（选定参考的ID名称）
# reference只有一条sequence
def ref_genome_size(file_path, reference_name):
    with open(file_path, 'r') as CC:
        # 读取第一行并去掉换行符
        current_line = CC.readline().strip()
        # 当current_line不是空行时
        while current_line:
            # 找到参考基因组，开始读取序列
            if current_line.startswith('>' + reference_name):
                # 重置sequence
                sequence = ''
                # 读取下一行并去掉换行符（序列的第一行）
                next_line = CC.readline().strip()
                while next_line and not next_line.startswith('>'):
                    # 累加序列
                    sequence += next_line
                    # 读取下一行并去掉换行符
                    next_line = CC.readline().strip()
                # 返回序列长度
                return len(sequence)
            # 如果当前行不是参考基因组，则读取下一行
            else:
                current_line = CC.readline().strip()
        # 如果文件遍历完毕仍未找到参考基因组，则返回0
        return 0

# 定义函数seqfromfasta_multi_genome，从包含多个基因组的fasta文件中，提取其中一个基因组的全部序列
# seqfromfasta_multi_genome（多基因组文件，想要提取的基因组id）
def seqfromfasta_multi_genome(multi_genome, split_genome_id, output_file):
    # 读取包含多个基因组序列的fasta文件
    DD = open(multi_genome, 'r')
    # 打开输出文件：用于写入需要提取的基因组序列
    EE = open(output_file, 'a+')
    # 按行读取
    DD_line = DD.readlines()
    # 遍历列表的每一个元素，即每一行
    for A in range(0, len(DD_line), 1):
        # 如果>和要提取的基因组名称在AA_line[A]中
        if '>' + split_genome_id in DD_line[A]:
            # 先写入fasta文件的首部
            EE.write(DD_line[A])
            # 再遍历后续行
            for B in range(A+1, len(DD_line), 1):
                # 如果出现下一个'>'，就退出循环
                if '>' in DD_line[B]:
                    break
                # 没有出现'>'或文件没有结束，将读取的序列写入输出文件
                else:
                    DD_content = DD_line[B].split('\n')[0]
                    EE.write(DD_content)
            EE.write('\n')
    DD.close()
    EE.close()

# 定义函数seqfromfasta_single_genome，从单个基因组的fasta格式文件中截取出某段序列。
# 定义函数seqfromfasta_single_genome(文件名称，起始位置，终止位置，输出文件名称）
def seqfromfasta_single_genome(single_genome, begin_position, end_position, output_file):
    begin_position = int(begin_position)
    end_position = int(end_position)
    # 读取single_genome文件
    FF = open(single_genome, 'r')
    # 打开输出文件
    GG = open(output_file, 'w')
    # 按行读取
    FF_line = FF.readlines()
    # 先写入首行
    GG.write(FF_line[0].strip('\n') + ' ' + str(begin_position) + '-' + str(end_position) + '\n')
    # 定义一个空列表
    sequence = []
    # 遍历CC_line的每一个元素，即single_genome的每一行
    for C in range(1, len(FF_line), 1):
        # 去除行末的换行符
        FF_content = FF_line[C].split('\n')[0]
        # 将FF_content作为元素，追加到sequence列表中
        sequence.append(FF_content)
    # 将sequence列表里的元素串接起来，没有间隔
    sequence = ''.join(sequence)
    # 通过索引的方式获得想要截取的序列，并写入输出文件
    GG.write(sequence[begin_position-1:end_position] + '\n')
    FF.close()
    GG.close()

def reverse_complement(sequence_file, output_file):
    # 定义碱基互补映射关系，包括大小写；不确定碱基N互补为自身
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n'
    }
    HH = open(sequence_file, 'r')
    II = open(output_file, 'w')
    sequence = ''
    HH_line = HH.readline()
    II.write(HH_line)
    HH_line = HH.readline()
    while HH_line:
        line_content = HH_line.split('\n')[0]
        sequence += line_content
        HH_line = HH.readline()
    # 反转序列并计算互补链
    # 反转序列，sequence[start:stop:step],step为1正向，为-1为负
    reversed_sequence = sequence[::-1]
    # 获得互补链累加赋值给reverse_complement_sequence
    reverse_complement_sequence = ''.join(complement.get(base, base) for base in reversed_sequence)
    II.write(reverse_complement_sequence + '\n')
    HH.close()
    II.close()

def main():
    # 创建解析器对象
    parser = argparse.ArgumentParser(description='The whole-mitogenome sequence alignment in fungi')
    # 添加必选参数
    parser.add_argument('--all_genome', required=True, help='a fasta file containing all genomes; Please modify the file name to end with .txt.')
    parser.add_argument('--reference_strainID', required=True, help="reference strain ID, such as 'DRR332962'")
    parser.add_argument('--block', type=int, required=True, help="The size of the smallest block, such as '200'")
    parser.add_argument('--cov', type=int, required=True, help='The smallest number of genomes aligned to the reference genome')
    # 解析命令行参数
    args = parser.parse_args()

    # 使用解析的参数
    all_genome = args.all_genome
    reference_strainID = args.reference_strainID
    block = args.block
    cov = args.cov

    # 后续代码可以继续使用 all_genome, reference_strainID, block, cov 变量
    print(f"输入参数: {all_genome}, {reference_strainID}, {block}, {cov}")

    # 定义参考基因组文件名称
    reference = '00'+'.'+reference_strainID+'.fa'
    # 获取参考基因组序列文件
    seqfromfasta_multi_genome(all_genome, reference_strainID, reference)
    # 使用参考基因组建库
    try:
        subprocess.run(['makeblastdb', '-in', reference, '-dbtype', 'nucl'])
    except subprocess.CalledProcessError as e:
        print(f"reference 建库 命令执行失败: {e}")
    # 将所有基因组与参考基因组两两比对
    try:
        subprocess.run(['blastn', '-db', reference, '-query', all_genome, '-out', '01.blastn.m6', '-evalue', '1e-5',
                        '-outfmt', '6', '-num_threads', '2'])
    except subprocess.CalledProcessError as e:
        print(f"blastn vs reference 命令执行失败: {e}")

    # 获取参考基因组的序列全长
    size = ref_genome_size(all_genome, reference_strainID)
    print('reference genome sequence length is ' + str(size))

    # 获取所有基因组的ID文件
    OA = open(all_genome, 'r')
    OB = open('00.genomeID_list', 'w')
    genome_file = OA.readline()
    while genome_file:
        if '>' in genome_file:
            genome_id = genome_file.split('\n')[0].split('>')[1]
            OB.write(genome_id + '\n')
        genome_file = OA.readline()
    OB.close()
    OA.close()

    # 记录参考基因组所有核苷酸位置，并设初始覆盖度为1
    dict_ref_all_position = {}
    for key in range(0, size, 1):
        dict_ref_all_position[key] = 1

    # 记录运行进度，输出便于知晓运行进展
    print('reading 01.blastn.m6', '......')

    # 统计所有基因组与参考基因组比对的结果，一个基因组比对参考基因组，参考基因组每个位置有且只记录覆盖1次
    OC = open('00.genomeID_list', 'r')
    OC_lines = OC.readlines()
    for ID in range(0, len(OC_lines), 1):
        processed_positions = {}
        OC_ID = OC_lines[ID].split('\n')[0]
        OD = open('01.blastn.m6', 'r')
        m6_list = OD.readlines()
        for i in range(0, len(m6_list), 1):
            m6_file = m6_list[i].split('\t')
            # 排除参考基因组自比
            if m6_file[0] != m6_file[1] and m6_file[0] == OC_ID:
                begin = m6_file[8]
                end = m6_file[9]
                if int(begin) > int(end):
                    begin = m6_file[9]
                    end = m6_file[8]
                for j in range(int(begin)-1, int(end), 1):
                    if processed_positions.get(j, 'no_value') == 'no_value':
                        dict_ref_all_position[j] = dict_ref_all_position[j] + 1
                        processed_positions[j] = 'already_recorded'
        OD.close()
    OC.close()

    # 这段代码的主要目的是根据前面统计的参考基因组碱基覆盖情况（存储在字典 dict 中），**********************************************
    # 找出覆盖度大于等于指定值 cov 且长度不小于指定值 block 的连续区域（即 block），**********************************************
    # 并将这些区域的起始和终止位置记录到文件 02.all_block 中。*****************************************************************
    OE = open('02.all_block', 'w')
    pos = int(-1)
    # 遍历字典中的键：遍历参考基因组的每一个基因
    for key in range(0, size, 1):
        # 如果键值大于所需参数cov的值，pos等于-1
        if dict_ref_all_position[key] >= cov and pos == -1:
            # pos等于键：pos等于基因所在位置（即一段block的起始位置）
            pos = key
        # 如果键值小于cov，pos不等于-1，且读段长度不少于cov，记录下参考基因组名称以及读段的起始终止位置
        elif dict_ref_all_position[key] < cov and pos != -1 and (key-1-pos+1) >= block:
            # 记录起始位置
            block_begin = pos+1
            # 记录终止位置
            block_end = key+1
            # 把这一个block写入文件02.all_block
            OE.write(reference_strainID+'\t'+str(block_begin)+'\t'+str(block_end)+'\n')
            # 更新pos的值重新为-1，寻找下一个block
            pos = int(-1)
        # 如果键值小于cov，pos不等于-1，且读段长度少于cov，舍弃直接重新寻找下一个block
        elif dict_ref_all_position[key] < cov and pos != -1 and (key-1-pos+1) < block:
            pos = int(-1)
    OE.close()

    # 对 02.all_block 文件中的 block 信息进行处理，当 block 的长度超过指定的上限（split_block）时，
    # 将其拆分成不超过该上限长度的子 block，并将拆分后的结果写入 03.all_block_split 文件。
    # block长度过长，以设定长度切分，如：999
    split_block = 999
    OF = open('02.all_block', 'r')
    OG = open('03.all_block_split', 'w')
    all_block_content = OF.readlines()
    for line in all_block_content:
        all_block_file = line.strip().split('\t')
        seq_id = str(all_block_file[0])
        seq_start = int(all_block_file[1])
        seq_end = int(all_block_file[2])
        while seq_start < seq_end:
            if seq_end - seq_start <= split_block:
                OG.write(seq_id + '\t' + str(seq_start) + '\t' + str(seq_end) + '\n')
                # 如果block长度小于或等于10000，退出检查下一个block
                break
            else:
                # 确保seq_mid不超过seq_end
                seq_mid = min(seq_start + split_block, seq_end)
                OG.write(seq_id + '\t' + str(seq_start) + '\t' + str(seq_mid) + '\n')
                # 移动到下一个块的开始位置
                seq_start = seq_mid + 1
    OF.close()
    OG.close()

    print('reading', all_genome, '......')

    # ******************************************创建一个字典：以基因ID为键，基因序列为键值*************************************
    genome_seq = {}
    one_genome_id = ''
    one_genome_seq = ''
    # 打开全部菌株的基因组fasta文件
    OH = open(all_genome, 'r')
    # 按行读取，获得一个列表文件
    all_genome_content = OH.readlines()
    # 遍历列表中的每一个元素，即文件的每一行；
    for k in range(0, len(all_genome_content), 1):
        # 如果存在'>'，且键值为空；就把这一行的基因ID存到键中
        if '>' in all_genome_content[k] and one_genome_seq == '':
            one_genome_id = all_genome_content[k].split(' ')[0]
            one_genome_id = one_genome_id.split('>')[1]
            one_genome_id = one_genome_id.split('\n')[0]
        # 如果存在'>'，且键值不为空(上一个键值刚结束)；重新将键值定义为空；再将这一行的基因ID存到键中
        elif '>' in all_genome_content[k] and one_genome_seq != '':
            genome_seq[one_genome_id] = one_genome_seq
            one_genome_seq = ''
            one_genome_id = all_genome_content[k].split(' ')[0]
            one_genome_id = one_genome_id.split('>')[1]
            one_genome_id = one_genome_id.split('\n')[0]
        # 如果这一行或者这个元素不存在'>'(即为序列行)，把序列累加在一起，存为键值
        else:
            one_genome_seq=str(one_genome_seq)+str(all_genome_content[k])
    # 如果键(one_genome_id)不在字典中，那么该键的值即位序列累加
    check_genome_id = genome_seq.get(one_genome_id, 'no_id')
    if check_genome_id == 'no_id':
        genome_seq[one_genome_id] = one_genome_seq
    OH.close()

    # 测试输出
    TT = open('Test_output.txt', 'w')

    # *************************共有的block再分别与所有基因组进行序列比对，获取所有基因组比对的上的区域*****************************
    print('reading', '03.all_block_split', '......')
    # 读取多基因组共有的block文件
    OI = open('03.all_block_split', 'r')
    # 获取一个按行读取的列表
    Input_file_content = OI.readlines()
    # 遍历列表中的元素：读取列表的每一行
    for i in range(0, len(Input_file_content), 1):
        # 剪切掉换行符
        Input_file_content[i] = Input_file_content[i].split('\n')[0]

        TT.write(Input_file_content[i])
        TT.write('\n')

        # 剪切'Tab'
        file = Input_file_content[i].split('\t')
        # 多基因共有block与参考基因组成功比对的起始位置
        subject_begin = file[1]
        # 多基因共有block与参考基因组成功比对的终止位置
        subject_end = file[2]

        # 定义参考 block 的文件名称
        ref_block_file = reference_strainID + '.' + str(subject_begin) + '-' + str(subject_end) + '.fa'
        # 获取参考 block 的序列文件
        seqfromfasta_single_genome(reference, subject_begin, subject_end, ref_block_file)
        # 参考 block 序列文件建库
        try:
            subprocess.run(['makeblastdb', '-in', ref_block_file, '-dbtype', 'nucl'])
        except subprocess.CalledProcessError as e:
            print(f"ref_block_file 建库 命令执行失败: {e}")
        # 所有基因组与参考 block 序列比对
        try:
            subprocess.run(['blastn', '-db', ref_block_file, '-query', all_genome, '-out', '04.blastn.m6', '-evalue', '1e-5',
                            '-outfmt', '6', '-num_threads', '2'])
        except subprocess.CalledProcessError as e:
            print(f"blastn vs ref_block_file 命令执行失败: {e}")

        # 定义单个 block 比对各基因组的序列输出文件名称
        one_block_ref_block_file = ('one' + '-' + 'block' + '-' + reference_strainID+'.' + str(subject_begin) + '-' +
                                    str(subject_end) + '.fa')

        # 判断是否已经存在该文件，如果存在，删去这个文件
        if os.path.exists(one_block_ref_block_file):
            try:
                subprocess.run(['rm', one_block_ref_block_file])
            except subprocess.CalledProcessError as e:
                print(f"rm one_block_ref_block_file 命令执行失败: {e}")
        # 参考 block 序列直接写入文件
        copy_file(ref_block_file, one_block_ref_block_file)

        TT.write('\n')

        print('reading', '04.blastn.m6', '......')

        # 对单个 block 比对后的文件进行筛选，选取identity>80,E值小于1e-50,
        OJ = open('04.blastn.m6', 'r')
        OK = open('05.select_blastn.m6', 'w')
        Input_content = OJ.readlines()
        for a in range(0, len(Input_content), 1):
            TT.write('****' + Input_content[a])
            data_line = Input_content[a].strip().split('\t')
            if data_line[0] != data_line[1]:
                if float(data_line[2]) > 80.0:
                    OK.write(Input_content[a])
        OJ.close()
        OK.close()

        TT.write('\n')

        # *****************读取筛选后的结果（05.select_identity_len_E_value.blastn.m6）按照第七列进行排序*********************
        # 如果该段block在一个基因组有多段比对上了，06.order_file就会得到一个排序后的文件06.order_file
        # 根据第七列排序（第七列是query的起始位置），就会将该段block在一个query基因组中的多个比对进行先后排序
        OL = open('05.select_blastn.m6', 'r')
        lines = OL.readlines()
        # 按照第七列进行排序，从小到大
        lines.sort(key=lambda line: int(line.split('\t')[6]))
        # lines.sort(key=lambda line: int(line.split('\t')[6]), reverse=True)
        for line in lines:
            line = line.split('\n')[0]
            TT.write(line + '*****' + '\n')
        # 排序后的文件写进文件06.order_file
        OM = open('06.order_file', 'w')
        OM.writelines(lines)
        OM.close()
        OL.close()

        TT.write('\n')

        #****************读取排序后的文件06.order_file，如果query基因组与该段block只比对上一段，这一段会直接保留*******************
        #****************如果query基因组与该段block有多段，第一段保留，如果第二段与第一段距离小于2000，保留第二段，以此类推***********
        #****************如果第二段与第一段距离过远，第二段舍弃，后面的第三段等也舍弃*********************************************
        line_seq = {}
        line_col7 = {}
        line_col8 = {}
        ON = open('06.order_file', 'r')
        line_content = ON.readlines()
        for i in range(0, len(line_content), 1):
            TT.write('******' + line_content[i])
            line_file = line_content[i].strip().split('\t')
            line_q_id = line_file[0]
            q_line_start = int(line_file[6])
            q_line_end = int(line_file[7])

            if line_seq.get(line_q_id, 'no value') == 'no value':
                # 如果字典中不存在该query基因组，直接以该基因组id为键，这一行的内容为值
                line_seq[line_q_id] = line_content[i]
                # 记录query基因组比对上的起始位置
                line_col7[line_q_id] = q_line_start
                # 记录query基因组比对上的终止位置
                line_col8[line_q_id] = q_line_end
            else:
                # 如果字典中已经记录了该query基因组，更新起始位置
                line_col7[line_q_id] = q_line_start
                # 如果更新后的起始位置大于等于此前记录的终止位置，并且差值小于设置范围（如：2000）
                if int(q_line_start) >= int(line_col8[line_q_id]) and abs(int(q_line_start) - int(line_col8[line_q_id])) <= 2000:
                    # 该query基因组比对上block的结果都保留
                    line_seq[line_q_id] = line_seq[line_q_id] + line_content[i]
                    # 更新query基因组比对上的起始位置
                    line_col7[line_q_id] = q_line_start
                    # 更新query基因组比对上的终止位置
                    line_col8[line_q_id] = q_line_end
                # 如果更新后的起始位置大于等于此前记录的终止位置，并且差值大于设置范围（如：2000）
                elif int(q_line_start) >= int(line_col8[line_q_id]) and abs(int(q_line_start) - int(line_col8[line_q_id])) >= 2000:
                    # query基因组比对上block的结果保持字典中记录的不变
                    line_seq[line_q_id] = line_seq[line_q_id]
                    # query基因组比对上的起始位置保持字典中里记录的不变
                    line_col7[line_q_id] = line_col7[line_q_id]
                    # query基因组比对上的终止位置保持字典中里记录的不变
                    line_col8[line_q_id] = line_col8[line_q_id]
        ON.close()

        TT.write('\n')

        # 把筛选后的结果写入06.order_file_result，作为排序后的最终结果
        OP = open('06.order_file_result', 'w')
        for id in line_seq:
            OP.write(line_seq[id])
            TT.write('*******' + line_seq[id])
        OP.close()

        TT.write('\n')

        # *********************读取06.order_file_result，获得新的block与各基因组比对的结果07.blastn.m6.result****************
        query_column1 = {}
        query_column1_start = {}
        query_column1_end = {}
        sub_start = {}
        sub_end = {}
        query_align_len = {}

        OQ = open('06.order_file_result', 'r')
        OR = open('07.blastn.m6.result', 'w')
        select_content = OQ.readlines()
        for i in range(0, len(select_content), 1):
            # 读取每行，获得列表
            select_file = select_content[i].strip().split('\t')
            # query基因组的名称
            ID_name = str(select_file[0])
            # query基因组比对上的序列长度
            select_sequence = int(select_file[3])
            # query基因组比对上block的起始位置
            select_q_start = int(select_file[6])
            # query基因组比对上block的终止位置
            select_q_end = int(select_file[7])
            # 参考block比对上的起始位置
            select_sub_start = int(select_file[8])
            # 参考block比对上的终止位置
            select_sub_end = int(select_file[9])

            # 如果query基因组不在query_column1中，即不存在该键
            if query_column1.get(ID_name, 'no value') == 'no value':
                # 在字典query_column1中记录该键，键值为空
                query_column1[ID_name] = ''
                # 在query_align_len中，记录query基因组比对上的序列长度
                query_align_len[ID_name] = select_sequence
                # 在query_column1_start中，记录query基因组比对上的起始位置
                query_column1_start[ID_name] = select_q_start
                # 在query_column1_end中，记录query基因组比对上的终止位置
                query_column1_end[ID_name] = select_q_end
                # 在sub_start中记录与参考block比对上的起始位置
                sub_start[ID_name] = select_sub_start
                # 在sub_end中记录与参考block比对上的终止位置
                sub_end[ID_name] = select_sub_end
            # 如果query基因组在query_column1中，即存在该键
            else:
                # 比对上的序列长度=字典中记录的+新一段的比对长度
                query_align_len[ID_name] = int(query_align_len[ID_name]) + int(select_sequence)
                # query起始位置=字典记录的起始位置与新一段的起始位置的最小值
                query_column1_start[ID_name] = min(int(query_column1_start[ID_name]), int(select_q_start))
                # query终止位置=字典记录的终止位置与新一段的终止位置的最大值
                query_column1_end[ID_name] = max(int(query_column1_end[ID_name]), int(select_q_end))
                if int(select_sub_start) < int(select_sub_end):
                    # 参考block的起始位置为所有比对段的起始位置的最小值
                    sub_start[ID_name] = min(int(sub_start[ID_name]), int(select_sub_start))
                    # 参考block的终止位置为所有比对段的终止位置的最大值
                    sub_end[ID_name] = max(int(sub_end[ID_name]), int(select_sub_end))
                else:
                    sub_start[ID_name] = max(int(sub_start[ID_name]), int(select_sub_start))
                    sub_end[ID_name] = min(int(sub_end[ID_name]), int(select_sub_end))

        # 记录下基因组比对的上block的结果，能比对上的query基因组都在query_column1键中
        for id in query_column1:
            # 自定义输出格式，部分值是不必需的，自定义数值
            OR.write(str(id) + '\t' + reference_strainID + '\t' + '90.317' + '\t' + str(query_align_len[id]) + '\t' + '1' +
                 '\t' + '1' + '\t' + str(query_column1_start[id]) + '\t' + str(query_column1_end[id]) + '\t' +
                 str(sub_start[id]) + '\t' + str(sub_end[id]) + '\t' + '0.0' + '\t' + '0' + '\n')
            TT.write(str(id) + '\t' + reference_strainID + '\t' + '90.317' + '\t' + str(query_align_len[id]) + '\t' + '1' +
                 '\t' + '1' + '\t' + str(query_column1_start[id]) + '\t' + str(query_column1_end[id]) + '\t' +
                 str(sub_start[id]) + '\t' + str(sub_end[id]) + '\t' + '0.0' + '\t' + '0' + '\n')
        OQ.close()
        OR.close()

        TT.write('\n')

        # ***********************读取07.blastn.m6.result，获取各基因组比对上参考block的读段，并用muscle比对************************
        print('reading', '07.blastn.m6.result', '...')
        # 读取07.blastn.m6.result文件
        OS = open('07.blastn.m6.result', 'r')
        # 按行读取，获得一个列表
        block_content = OS.readlines()
        # 遍历列表的每一个元素：读取文件的每一行
        for j in range(0, len(block_content), 1):
            # 剪切掉Tab
            block_content_split = block_content[j].split('\t')
            # 文件第一列为基因组名称
            query_id = block_content_split[0]
            # 基因组与参考block读段比对上的起始位置
            query_begin = block_content_split[6]
            # 基因组与参考block读段比对上的终止位置
            query_end = block_content_split[7]
            # 避免相同基因组进行比对（但好像不会和自己比较，可能多余了）
            if block_content_split[0] != block_content_split[1]:
                # 剪切出比对得上的序列：得到block_align文件（fasta格式）
                # 先剪切出特定的基因组
                seqfromfasta_multi_genome(all_genome, block_content_split[0], block_content_split[0] + '.fa')
                # 定义block_align文件
                block_align = block_content_split[0] + '.' + str(query_begin) + '-' + str(query_end) + '.fa'
                # 再进行特定位置序列的剪切
                seqfromfasta_single_genome(block_content_split[0] + '.fa', query_begin, query_end, block_align)
                # 判断与参考基因组是否正向比对，如果不是就取反向互补序列
                if int(block_content_split[8]) < int(block_content_split[9]):
                    copy_file(block_align, one_block_ref_block_file)
                else:
                    # 定义反向互补链文件
                    rev_file = block_content_split[0] + '.' + query_begin + '-' + query_end + '_rev.fa'
                    # 取反向互补连
                    reverse_complement(block_align, rev_file)
                    # 追加反向互补序列文件到one_block_ref_block_file
                    copy_file(block_content_split[0] + '.' + query_begin + '-' + query_end + '_rev.fa',
                              one_block_ref_block_file)
        OS.close()

        # 对one_block_ref_block_file进行muscle比对
        # 定义生成的muscle文件名称
        one_block_ref_block_file_muscle = ('one'+'-'+'block'+'-'+reference_strainID+'.'+str(subject_begin)+'-' +
                                           str(subject_end)+'.muscle')
        # 运行muscle比对
        try:
            subprocess.run(['muscle', '-align', one_block_ref_block_file, '-output', one_block_ref_block_file_muscle],
                           check=True)
        # muscle命令操作异常处理
        except subprocess.CalledProcessError as e:
            print(f"muscle 命令执行失败: {e}")

        # **************由于部分block并不是所有基因组都具备，因此补充没被记录的基因组，并用'-----'表示序列长度**********************
        # **************如果有12个基因组，选择cov为7,那么会筛选出有些block比对上的参考基因组会没有序列，使用'-'代替******************
        # 下面就是记录所有比对上的基因组ID，以及muscle比对后的序列长度。再与所有基因组比较，发现哪些基因组没有比对上，就补充'-'***********
        # 创建一个字典：ID为键，序列为值
        block_seq_id = {}
        # 定义一个空ID名称
        block_one_seq_id = ''
        print('reading', one_block_ref_block_file_muscle, '...')
        # 读取one_block_ref_block_file_muscle文件
        OT = open(one_block_ref_block_file_muscle, 'r')
        # 按行读取文件
        muscle_content = OT.readlines()
        # 遍历文件的每一行
        for z in range(0, len(muscle_content), 1):
            # 如果'>'存在，记录基因组ID
            if '>' in muscle_content[z]:
                block_one_seq_id = muscle_content[z].split(' ')[0]
                block_one_seq_id = block_one_seq_id.split('>')[1]
            else:
                # 如果字典中不存在这个键(基因组ID),记录该键的值为0
                if block_seq_id.get(block_one_seq_id, 'no_exist') == 'no_exist':
                    block_seq_id[block_one_seq_id] = 0
                # 切掉该行的换行符(因为要记录序列的长度)
                muscle_content[z] = muscle_content[z].split('\n')[0]
                # 记录序列长度：序列长度累加
                block_seq_id[block_one_seq_id] = block_seq_id[block_one_seq_id]+len(muscle_content[z])
        OT.close()

        # 遍历genome_seq中的键，该字典具有所有基因组的ID名称
        for id in genome_seq:
            # 如果one_block_ref_block_file_muscle文件没有该基因组ID
            if block_seq_id.get(id, 'no_id') == 'no_id':
                # 在one_block_ref_block_file_muscle文件后追加
                OU = open(one_block_ref_block_file_muscle, 'a+')
                # 追加补充的'>基因组ID名称'以及换行符
                OU.write('>' + id + '\n')
                b_length = 0
                # 遍历字典block_seq_id的键
                for b_id in block_seq_id:
                    # 记录键值：序列长度
                    b_length = block_seq_id[b_id]
                    break
                # 追加同样字符长度的'-'
                OU.write(b_length * '-' + '\n')
                OU.close()
    OI.close()

    # 定义要搜索的目录（当前目录）
    directory = os.getcwd()
    # 列出当前目录下的所有文件
    for filename in os.listdir(directory):
        # 检查文件名是否包含.fa
        if '.fa' in filename:
            # 构建完整的文件路径
            filepath = os.path.join(directory, filename)
            # 检查是否是一个文件而不是目录
            if os.path.isfile(filepath):
                try:
                    os.remove(filepath)
                except OSError as e:
                    print('delete error')

if __name__ == "__main__":
    main()


