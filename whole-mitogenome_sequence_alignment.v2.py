# 引入模块，可以在python系统下运行linux命令
import os
import subprocess
import argparse
from collections import defaultdict

# 定义函数：把一个文件的内容追加到另一个文件中
def copy_file (file_1, file_2):
    # 读取文件1
    OE = open(file_1, 'r')
    # 按行读取
    file_1_content = OE.readlines()
    # 打开文件2，追加
    OF = open(file_2, 'a+')
    # 遍历列表
    for x in range(0, len(file_1_content), 1):
        # 将文件1的每行追加到文件2中
        OF.write(file_1_content[x])
    OE.close()
    OF.close()

#定义函数：读取FASTA格式文件获得特定ID的序列长度
#定义函数ref_genome_size，所需参数file_path（文件路径）；reference_name（选定参考的ID名称）
#reference只有一条sequence
def ref_genome_size(file_path, reference_name):
    with open(file_path, 'r') as file:
        #读取第一行并去掉换行符
        current_line = file.readline().strip()
        #当current_line不是空行时
        while current_line:
            #找到参考基因组，开始读取序列
            if current_line.startswith('>' + reference_name):
                #重置sequence
                sequence = ''
                #读取下一行并去掉换行符（序列的第一行）
                next_line = file.readline().strip()
                while next_line and not next_line.startswith('>'):
                    # 累加序列
                    sequence += next_line
                    # 读取下一行并去掉换行符
                    next_line = file.readline().strip()
                #返回序列长度
                return len(sequence)
            #如果当前行不是参考基因组，则读取下一行
            else:
                current_line = file.readline().strip()
        #如果文件遍历完毕仍未找到参考基因组，则返回0
        return 0

#定义函数seqfromfasta_multi_genome，从包含多个基因组的fasta文件中，提取其中一个基因组的全部序列
#seqfromfasta_multi_genome（多基因组文件，想要提取的基因组id）
def seqfromfasta_multi_genome(multi_genome, split_genome_id, output_file):
    #读取包含多个基因组序列的fasta文件
    AA = open(multi_genome, 'r')
    #打开输出文件：用于写入需要提取的基因组序列
    BB = open(output_file, 'a+')
    #按行读取
    AA_line = AA.readlines()
    #遍历列表的每一个元素，即每一行
    for A in range(0, len(AA_line), 1):
        #如果>和要提取的基因组名称在AA_line[A]中
        if '>' + split_genome_id in AA_line[A]:
            #先写入fasta文件的首部
            BB.write(AA_line[A])
            #再遍历后续行
            for B in range(A+1, len(AA_line), 1):
                #如果出现下一个'>'，就退出循环
                if '>' in AA_line[B]:
                    break
                #没有出现'>'或文件没有结束，将读取的序列写入输出文件
                else:
                    AA_content = AA_line[B].split('\n')[0]
                    BB.write(AA_content)
            BB.write('\n')
    AA.close()
    BB.close()

#定义函数seqfromfasta_single_genome，从单个基因组的fasta格式文件中截取出某段序列。
#定义函数seqfromfasta_single_genome(文件名称，起始位置，终止位置，输出文件名称）
def seqfromfasta_single_genome(single_genome, begin_position, end_position, output_file):
    begin_position = int(begin_position)
    end_position = int(end_position)
    #读取single_genome文件
    CC = open(single_genome, 'r')
    #打开输出文件
    EE = open(output_file, 'w')
    #按行读取
    CC_line = CC.readlines()
    #先写入首行
    EE.write(CC_line[0].strip('\n') + ' ' + str(begin_position) + '-' + str(end_position) + '\n')
    #定义一个空列表
    sequence = []
    #遍历CC_line的每一个元素，即single_genome的每一行
    for C in range(1, len(CC_line), 1):
        #去除行末的换行符
        CC_content = CC_line[C].split('\n')[0]
        #将CC_content作为元素，追加到sequence列表中
        sequence.append(CC_content)
    #将sequence列表里的元素串接起来，没有间隔
    sequence = ''.join(sequence)
    #通过索引的方式获得想要截取的序列，并写入输出文件
    EE.write(sequence[begin_position-1:end_position] + '\n')
    CC.close()
    EE.close()

def seqfromfasta_single_genome_just_sequence(single_genome, begin_position, end_position, output_file):
    begin_position = int(begin_position)
    end_position = int(end_position)
    # 读取single_genome文件
    CC = open(single_genome, 'r')
    # 打开输出文件
    EE = open(output_file, 'w')
    # 按行读取
    CC_line = CC.readlines()
    # 定义一个空列表
    sequence = []
    # 遍历CC_line的每一个元素，即single_genome的每一行
    for C in range(1, len(CC_line), 1):
        # 去除行末的换行符
        CC_content = CC_line[C].split('\n')[0]
        # 将CC_content作为元素，追加到sequence列表中
        sequence.append(CC_content)
    # 将sequence列表里的元素串接起来，没有间隔
    sequence = ''.join(sequence)
    # 通过索引的方式获得想要截取的序列，并写入输出文件
    EE.write(sequence[begin_position - 1:end_position])
    CC.close()
    EE.close()

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

    #根据选择的参考基因组建库
    reference = '00'+'.'+reference_strainID+'.fa'
    seqfromfasta_multi_genome(all_genome, reference_strainID, reference)
    subprocess.run(['makeblastdb', '-in', reference, '-dbtype', 'nucl'])
    subprocess.run(['blastn', '-db', reference, '-query', all_genome, '-out', '01.blastn.m6', '-evalue', '1e-5', '-outfmt',
                    '6', '-num_threads', '2'])
    
    #根据定义的函数记录参考基因组序列长度
    size = ref_genome_size(all_genome, reference_strainID)
    print('reference genome sequence length is ' + str(size))
    
    #记录所有基因组的ID
    OY = open(all_genome, 'r')
    OX = open('00.genomeID_list', 'w')
    genome_file = OY.readline()
    while genome_file:
        if '>' in genome_file:
            genome_id = genome_file.split('\n')[0].split('>')[1]
            OX.write(genome_id + '\n')
        genome_file = OY.readline()
    OX.close()
    OY.close()
    
    #根据比对结果，记录参考基因组每个核苷酸位置的覆盖度
    dict = {}
    for key in range(0, size, 1):
        dict[key] = 1
    
    print('reading 01.blastn.m6_new', '01.blastn.m6_new', '...')
    DL = open('00.genomeID_list', 'r')
    DL_lines = DL.readlines()
    for ID in range(0, len(DL_lines), 1):
        processed_positions = {}
        DL_ID = DL_lines[ID].split('\n')[0]
        F = open('01.blastn.m6', 'r')
        m6_list = F.readlines()
        for i in range(0, len(m6_list), 1):
            m6_file = m6_list[i].split('\t')
            if m6_file[0] != m6_file[1] and m6_file[0] == DL_ID:
                begin = m6_file[8]
                end = m6_file[9]
                if int(begin) > int(end):
                    begin = m6_file[9]
                    end = m6_file[8]
                for j in range(int(begin)-1, int(end), 1):
                    if processed_positions.get(j, 'no_value') == 'no_value':
                        dict[j] = dict[j] + 1
                        processed_positions[j] = 'already_recorded'
        F.close()
    DL.close()

    print(dict.values())

    #这段代码的主要目的是根据前面统计的参考基因组碱基覆盖情况（存储在字典 dict 中），***************************************************
    #找出覆盖度大于等于指定值 cov 且长度不小于指定值 block 的连续区域（即 block），***************************************************
    #并将这些区域的起始和终止位置记录到文件 02.all_block 中。**********************************************************************
    OH = open('02.all_block', 'w')
    pos = int(-1)
    #遍历字典中的键：遍历参考基因组的每一个基因
    for key in range(0, size, 1):
        #如果键值大于所需参数cov的值，pos等于-1
        if dict[key] >= cov and pos == -1:
            #pos等于键：pos等于基因所在位置（即一段block的起始位置）
            pos = key
        #如果键值小于cov，pos不等于-1，且读段长度不少于cov，记录下参考基因组名称以及读段的起始终止位置
        elif dict[key] < cov and pos != -1 and (key-1-pos+1) >= block:
            #记录起始位置
            block_begin = pos+1
            #记录终止位置
            block_end = key+1
            #把这一个block写入文件02.all_block
            OH.write(reference_strainID+'\t'+str(block_begin)+'\t'+str(block_end)+'\n')
            #更新pos的值重新为-1，寻找下一个block
            pos = int(-1)
        # 如果键值小于cov，pos不等于-1，且读段长度少于cov，舍弃直接重新寻找下一个block
        elif dict[key] < cov and pos != -1 and (key-1-pos+1) < block:
            pos = int(-1)
    OH.close()
    
    #对 02.all_block 文件中的 block 信息进行处理，当 block 的长度超过指定的上限（split_block）时，
    #将其拆分成不超过该上限长度的子 block，并将拆分后的结果写入 03.all_block_split 文件。
    #block长度过长，以设定长度切分，如：999
    split_block = 2999
    OP = open('02.all_block', 'r')
    OQ = open('03.all_block_split', 'w')
    all_block_content = OP.readlines()
    for line in all_block_content:
        all_block_file = line.strip().split('\t')
        seq_id = all_block_file[0]
        seq_start = int(all_block_file[1])
        seq_end = int(all_block_file[2])
        while seq_start < seq_end:
            if seq_end - seq_start <= split_block:
                OQ.write(seq_id + '\t' + str(seq_start) + '\t' + str(seq_end) + '\n')
                # 如果block长度小于或等于10000，退出检查下一个block
                break
            else:
                # 确保seq_mid不超过seq_end
                seq_mid = min(seq_start + split_block, seq_end)
                OQ.write(seq_id + '\t' + str(seq_start) + '\t' + str(seq_mid) + '\n')
                # 移动到下一个块的开始位置
                seq_start = seq_mid + 1
    OP.close()
    OQ.close()
    
    #******************************************创建一个字典：以基因ID为键，基因序列为键值*****************************************
    genome_seq = {}
    one_genome_id = ""
    one_genome_seq = ""
    print('reading', all_genome, '...')
    #打开全部菌株的基因组fasta文件
    OC = open(all_genome, 'r')
    #按行读取，获得一个列表文件
    all_genome_content = OC.readlines()
    #遍历列表中的每一个元素，即文件的每一行；
    for k in range(0, len(all_genome_content), 1):
        #如果存在'>'，且键值为空；就把这一行的基因ID存到键中
        if ">" in all_genome_content[k] and one_genome_seq == "":
             one_genome_id = all_genome_content[k].split(" ")[0]
             one_genome_id = one_genome_id.split('>')[1]
             one_genome_id = one_genome_id.split('\n')[0]
        #如果存在'>'，且键值不为空(上一个键值刚结束)；重新将键值定义为空；再将这一行的基因ID存到键中
        elif ">" in all_genome_content[k] and one_genome_seq != "":
             genome_seq[one_genome_id] = one_genome_seq
             one_genome_seq = ""
             one_genome_id = all_genome_content[k].split(" ")[0]
             one_genome_id = one_genome_id.split('>')[1]
             one_genome_id = one_genome_id.split('\n')[0]
        #如果这一行或者这个元素不存在'>'(即为序列行)，把序列累加在一起，存为键值
        else:
             one_genome_seq = str(one_genome_seq)+str(all_genome_content[k])
    #如果键(one_genome_id)不在字典中，那么该键的值即位序列累加
    check_genome_id = genome_seq.get(one_genome_id, 'no_id')
    if check_genome_id == 'no_id':
        genome_seq[one_genome_id] = one_genome_seq
    OC.close()
    
    #测试输出
    TT = open('Test_output.txt', 'w')
    
    #*************************共有的block再分别与所有基因组进行序列比对，获取所有基因组比对的上的区域**********************************
    print('reading', '03.all_block_split', '...')
    #读取多基因组共有的block文件
    OA = open('03.all_block_split', 'r')
    #获取一个按行读取的列表
    Input_file_content = OA.readlines()
    #遍历列表中的元素：读取列表的每一行
    for i in range(0, len(Input_file_content), 1):
        #剪切掉换行符
        Input_file_content[i] = Input_file_content[i].split('\n')[0]
    
        TT.write(Input_file_content[i])
        TT.write('\n')
    
        #剪切'Tab'
        file = Input_file_content[i].split('\t')
        #多基因共有block与参考基因组成功比对的起始位置
        subject_begin = file[1]
        #多基因共有block与参考基因组成功比对的终止位置
        subject_end = file[2]
        #记录共有block的长度(自己加的）
        block_shared_length = int(subject_end)-int(subject_begin)+1
        # print(block_shared_length)
        #剪切出block序列进行建库，并与多个基因组进行比对，获取一个block的所有基因组比对的合并文件
        ref_block_file = reference_strainID + '.' + str(subject_begin) + '-' + str(subject_end) + '.fa'
    
        seqfromfasta_single_genome(reference, subject_begin, subject_end, ref_block_file)
        #subprocess.run(['perl', seqfromfasta, reference, reference_strainID, ref_block_file, subject_begin, subject_end])
        subprocess.run(['makeblastdb', '-in', ref_block_file, '-dbtype', 'nucl'])
        subprocess.run(['blastn', '-db', ref_block_file, '-query', all_genome, '-out', '04.blastn.m6', '-evalue', '1e-5',
                        '-outfmt', '6', '-num_threads', '2'])
        one_block_ref_block_file = ('one' + '-' + 'block' + '-' + reference_strainID+'.' + str(subject_begin) + '-' +
                                    str(subject_end) + '.fa')
    
        #判断是否已经存在该文件，如果存在，删去这个文件
        if os.path.exists(one_block_ref_block_file):
            subprocess.run(['rm', one_block_ref_block_file])
        #使用定义函数将ref_block_file写入one_block_ref_block_file中:把参考的block序列写进one_block_ref_block_file
        copy_file(ref_block_file, one_block_ref_block_file)
    
        TT.write('\n')
    
        #**********对一段block比对后的文件（04.blastn.m6）进行筛选，选取长度大于block，identity>85,E值小于1e-50,*********************
        print('reading', '04.blastn.m6', '...')
        OM = open('04.blastn.m6', 'r')
        ON = open('05.select_identity_len_E_value.blastn.m6', 'w')
        Input_content = OM.readlines()
        for a in range(0, len(Input_content), 1):
            TT.write('****' + Input_content[a])
            data_line = Input_content[a].strip().split('\t')
            if data_line[0] != data_line[1]:
                #if (float(data_line[2]) > 70.0 and int(data_line[3]) > block_shared_length*0.8 and float(data_line[10])
                #        < 1.0e-50):
                #if int(data_line[3]) > block and float(data_line[10]) < 1.0e-50:
                #if float(data_line[2]) > 70.0 and int(data_line[3]) > block_shared_length*0.8:
                #做了多次测试80是较好的选择
                if float(data_line[2]) > 80.0 and float(data_line[10]) <= 1e-50:
                    ON.write(Input_content[a])
        OM.close()
        ON.close()
    
        TT.write('\n')
    
        #*****************读取筛选后的结果（05.select_identity_len_E_value.blastn.m6）*****************************************
        #记录每个基因组比对该段block的最大长度的那一段
    
        DA = open('05.select_identity_len_E_value.blastn.m6', 'r')
        lines = DA.readlines()
        for line in lines:
            line = line.split('\n')[0]
            TT.write(line + '*****' + '\n')
    
        TT.write('\n')
    
        #创建字典，id为键，记录最大的比对长度
        q_align_4_lie_max = {}
        #创建字典，以不同基因组id为键，‘最大的比对长度的那一行’为值
        q_align_4_lie_max_line = {}
    
        XX = open('00.genomeID_list', 'r')
        XX_id_list = XX.readlines()
        for i in range(0, len(XX_id_list), 1):
            XX_id = XX_id_list[i].split('\n')[0]
            YY = open('05.select_identity_len_E_value.blastn.m6', 'r')
            YY_line_list = YY.readlines()
            for t in range(0, len(YY_line_list), 1):
                YY_line = YY_line_list[t]
                YY_line_id = YY_line_list[t].split('\n')[0].split('\t')[0]
                refer_id = YY_line_list[t].split('\n')[0].split('\t')[1]
                lie_4_length = int(YY_line_list[t].split('\n')[0].split('\t')[3])
                if YY_line_id != refer_id and YY_line_id == XX_id:
                    if XX_id not in q_align_4_lie_max:
                        q_align_4_lie_max[XX_id] = lie_4_length
                        q_align_4_lie_max_line[XX_id] = YY_line
                    if XX_id in q_align_4_lie_max and lie_4_length > q_align_4_lie_max[XX_id]:
                        q_align_4_lie_max[XX_id] = lie_4_length
                        q_align_4_lie_max_line[XX_id] = YY_line
            YY.close()
        XX.close()
    
        #print(q_align_4_lie_max_line.values())
        #找到比对最长的那一段，在这一段上下浮动一定长度寻找block；如：2000
        GG = open('06.select_line_file', 'w')
        for id in q_align_4_lie_max_line:
            max_start = int(q_align_4_lie_max_line[id].split('\n')[0].split('\t')[6])
            max_end = int(q_align_4_lie_max_line[id].split('\n')[0].split('\t')[7])
            select_start = int(max_start) - 3000
            select_end = int(max_end + 3000)
            #考虑到最长的部分可能在首部，或者尾部。范围要重新划分
            #如果首部-2000的值为负值，则原范围分为两部分：（1，select_end）和（select_start2，ref_genome_size(all_genome, id)）
            if select_start <= 0:
                select_start1 = 1
                select_start2 = ref_genome_size(all_genome, id) + select_start
                ZZ = open('05.select_identity_len_E_value.blastn.m6', 'r')
                ZZ_line = ZZ.readlines()
                for i in range(0, len(ZZ_line), 1):
                    ZZ_content = ZZ_line[i]
                    ZZ_id = ZZ_line[i].split('\n')[0].split('\t')[0]
                    ZZ_q_start = int(ZZ_line[i].split('\n')[0].split('\t')[6])
                    ZZ_q_end = int(ZZ_line[i].split('\n')[0].split('\t')[7])
                    if id == ZZ_id and ZZ_q_start >= 1 and ZZ_q_end <= select_end:
                        GG.write(ZZ_content)
                        #print(ZZ_content)
                    if id == ZZ_id and ZZ_q_start >= select_start2 and ZZ_q_end <= ref_genome_size(all_genome, id):
                        GG.write(ZZ_content)
                        #print(ZZ_content)
                ZZ.close()
    
            #如果尾部+2000超出该基因组的总长的话，则原范围也要划分两部分：（select_start---ref_genome_size(all_genome, id)）和（1，select_end2）
            elif select_end > ref_genome_size(all_genome, id):
                select_end1 = ref_genome_size(all_genome, id)
                select_end2 = select_end - ref_genome_size(all_genome, id)
                ZZ = open('05.select_identity_len_E_value.blastn.m6', 'r')
                ZZ_line = ZZ.readlines()
                for i in range(0, len(ZZ_line), 1):
                    ZZ_content = ZZ_line[i]
                    ZZ_id = ZZ_line[i].split('\n')[0].split('\t')[0]
                    ZZ_q_start = int(ZZ_line[i].split('\n')[0].split('\t')[6])
                    ZZ_q_end = int(ZZ_line[i].split('\n')[0].split('\t')[7])
                    if id == ZZ_id and ZZ_q_start >= select_start and ZZ_q_end <= ref_genome_size(all_genome, id):
                        GG.write(ZZ_content)
                        print(ZZ_content)
                    if id == ZZ_id and ZZ_q_start >= 1 and ZZ_q_end <= select_end2:
                        GG.write(ZZ_content)
                        print(ZZ_content)
                ZZ.close()
    
            #这是最常见的结果就是；首部尾部+-2000均未出现首部小于0，尾部超出基因组长度。则比对上的block须在(select_start, select_end)
            elif select_start > 0 and select_end <= ref_genome_size(all_genome, id):
                #print(select_start, select_end)
                ZZ = open('05.select_identity_len_E_value.blastn.m6', 'r')
                ZZ_line = ZZ.readlines()
                for i in range(0, len(ZZ_line), 1):
                    ZZ_content = ZZ_line[i]
                    ZZ_id = ZZ_line[i].split('\n')[0].split('\t')[0]
                    ZZ_q_start = int(ZZ_line[i].split('\n')[0].split('\t')[6])
                    ZZ_q_end = int(ZZ_line[i].split('\n')[0].split('\t')[7])
                    if id == ZZ_id and ZZ_q_start >= select_start and ZZ_q_end <= select_end:
                        GG.write(ZZ_content)
                        #print(ZZ_content)
                ZZ.close()
        GG.close()
    
        WW = open('06.select_line_file', 'r')
        WW_lines = WW.readlines()
        for i in range(0, len(WW_lines), 1):
            TT.write('******' + WW_lines[i])
        TT.write('\n')
    
        #以第四列排序，从大到小
        WW_lines.sort(key=lambda line: int(line.split('\t')[3]), reverse=True)
        JJ = open('07.sort_4_file', 'w')
        JJ.writelines(WW_lines)
        JJ.close()
    
        # 以ID为键，对应的行为值；记录ID所有的行
        select = {}
        PP = open('00.genomeID_list', 'r')
        ID_line = PP.readlines()
        for i in range(0, len(ID_line), 1):
            ID_select = ID_line[i].split('\n')[0]
            DC = open('07.sort_4_file', 'r')
            DC_line = DC.readlines()
            for j in range(0, len(DC_line), 1):
                l_content = DC_line[j]
                l_ID = DC_line[j].split('\n')[0].split('\t')[0]
                if ID_select not in select and ID_select == l_ID:
                    select[ID_select] = [l_content]
                elif ID_select in select and ID_select == l_ID:
                    select[ID_select].append(l_content)
            DC.close()
        PP.close()
    
        QQ = open('07.sort_4_file_result', 'w')
        len_line = {}
        for id in select:
            #如果ID对应的只有一行
            if len(select[id]) <= 1:
                QQ.write(select[id][0])
            #如果ID对应有多行
            else:
                for g in range(0, len(select[id]), 1):
                    #读取每一行比对上的长度
                    length = select[id][g].split('\n')[0].split('\t')[3]
                    #提前转换为整数
                    length_int = int(length)
                    # 如果第一段的长度只在共有block长度向下浮动20bp或大于共有block长度；只要这一段其余的不要；退出循环
                    if len_line.get(id, 'no_value') == 'no_value' and length_int >= (block_shared_length - 10):
                        QQ.write(select[id][g])
                        break
    
                    # 如果第一段的长度短于共有block长度-20，这一行保留；记录在字典中
                    elif len_line.get(id, 'no_value') == 'no_value' and length_int <= (block_shared_length - 10):
                        # 字典中记录这一段的长度
                        len_line[id] = length_int
                        # 输出此行
                        QQ.write(select[id][g])
    
                    # 如果字典中存在,并且片段与字典中比对上的长度总和小于（block_shared_length - 20)；这一段也保留
                    elif len_line.get(id, 'no_value') != 'no_value' and (length_int + len_line[id]) <= (
                            block_shared_length - 10):
                        len_line[id] += length_int
                        QQ.write(select[id][g])
    
                    # 如果与字典中长度相加在参考block长度上下浮动20,这一段也要。
                    elif len_line.get(id, 'no_value') != 'no_value' and (block_shared_length - 10) <= (
                            length_int + len_line[id]) <= (block_shared_length + 10):
                        len_line[id] += length_int
                        QQ.write(select[id][g])
    
                    #
                    elif len_line.get(id, 'no_value') != 'no_value' and (length_int + len_line[id]) > (block_shared_length + 10):
                        QQ.write(select[id][g])
                        break
        QQ.close()
    
        #读取文件，先找与参考block正向比对的序列；再找反向比对的序列。分别写到两个文件中
        OK = open('07.sort_4_file_result', 'r')
        select_content = OK.readlines()
        #记录比对上参考block正向序列的行
        AB = open('08.lines_forward_align', 'w')
        #记录比对上参考block反向序列的行
        BC = open('08.lines_reverse_align', 'w')
        for i in range(0, len(select_content), 1):
            TT.write('#######' + select_content[i])
            ref_block_start = int(select_content[i].split('\n')[0].split('\t')[8])
            ref_block_end = int(select_content[i].split('\n')[0].split('\t')[9])
            if ref_block_start < ref_block_end:
                #正向比对的行写入文件
                AB.write(select_content[i])
            if ref_block_start > ref_block_end:
                BC.write(select_content[i])
        OK.close()
        AB.close()
        BC.close()
    
        # 对正向序列比对的行进行排序：第9列从小到大，然后依次读取，截取序列，追加到文件（）中
        CD = open('08.lines_forward_align', 'r')
        CD_lines = CD.readlines()
        CD_lines.sort(key=lambda line: int(line.split('\t')[8]))
        # print(CD_lines)
        DE = open('08.lines_forward_align_sort', 'w')
        DE.writelines(CD_lines)
        DE.close()
    
        print('\n')
    
        #对反向序列比对的行进行排序：第9列从大到小，然后依次读取，截取序列，追加到文件（）中
        EF = open('08.lines_reverse_align', 'r')
        EF_lines = EF.readlines()
        EF_lines.sort(key=lambda line: int(line.split('\t')[8]), reverse=True)
        # print(EF_lines)
        FG = open('08.lines_reverse_align_sort', 'w')
        FG.writelines(EF_lines)
        FG.close()
    
        #创建字典记录ID名称以及序列
        id_record = {}
        #读取正向比对排序后的文件，并剪切比对的序列
        HA = open(one_block_ref_block_file, 'a+')
        HI = open('00.genomeID_list')
        HI_list = HI.readlines()
        for n in range(0, len(HI_list), 1):
            id_name = HI_list[n].split('\n')[0]
            GH = open('08.lines_forward_align_sort', 'r')
            GH_lines = GH.readlines()
            for r in range(0, len(GH_lines), 1):
                q_ID = GH_lines[r].split('\n')[0].split('\t')[0]
                q_start = GH_lines[r].split('\n')[0].split('\t')[6]
                q_end = GH_lines[r].split('\n')[0].split('\t')[7]
                if id_name == q_ID and q_ID not in id_record:
                    HA.write('>' + q_ID + '\n')
                    HA.flush()
                    q_ID_file = q_ID + '.fa'
                    seqfromfasta_multi_genome(all_genome, q_ID, q_ID_file)
                    q_ID_split = q_ID + '.split.fa'
                    seqfromfasta_single_genome_just_sequence(q_ID_file, q_start, q_end, q_ID_split)
                    copy_file(q_ID_split, one_block_ref_block_file)
                    id_record[q_ID] = ''
                elif id_name == q_ID and q_ID in id_record:
                    q_ID_file = q_ID + '.fa'
                    #这个文件会覆盖
                    q_ID_split = q_ID + '.split.fa'
                    seqfromfasta_single_genome_just_sequence(q_ID_file, q_start, q_end, q_ID_split)
                    copy_file(q_ID_split, one_block_ref_block_file)
            HA.write('\n')
            HA.flush()
            GH.close()
        HI.close()
        HA.close()
    
        # 创建字典记录ID名称
        id_record_new = {}
        # 读取反向比对排序后的文件，并剪切比对的序列
        HC = open('00.genomeID_list')
        HC_list = HC.readlines()
        for m in range(0, len(HC_list), 1):
            id_m = HC_list[m].split('\n')[0]
            HB = open('reverse.file', 'w')
            HD = open('08.lines_reverse_align_sort', 'r')
            HD_lines = HD.readlines()
            for v in range(0, len(HD_lines), 1):
                qq_ID = HD_lines[v].split('\n')[0].split('\t')[0]
                qq_start = HD_lines[v].split('\n')[0].split('\t')[6]
                qq_end = HD_lines[v].split('\n')[0].split('\t')[7]
                #print(qq_start)
                #print(qq_end)
                if id_m == qq_ID and qq_ID not in id_record_new:
                    HB.write('>' + qq_ID + '\n')
                    HB.flush()
                    qq_ID_file = qq_ID + '.fa'
                    seqfromfasta_multi_genome(all_genome, qq_ID, qq_ID_file)
                    qq_ID_split = qq_ID + '.split.fa'
                    seqfromfasta_single_genome_just_sequence(qq_ID_file, qq_start, qq_end, qq_ID_split)
                    copy_file(qq_ID_split, 'reverse.file')
                    id_record_new[qq_ID] = ''
                elif id_m == qq_ID and qq_ID in id_record_new:
                    qq_ID_file = qq_ID + '.fa'
                    # 这个文件会覆盖
                    qq_ID_split = qq_ID + '.split.fa'
                    seqfromfasta_single_genome_just_sequence(qq_ID_file, qq_start, qq_end, qq_ID_split)
                    copy_file(qq_ID_split, 'reverse.file')
            HB.write('\n')
            HB.flush()
            HB.close()
            # subprocess.run([reverse_perl, 'reverse.file', 'reverse.file_rev'])
            reverse_complement('reverse.file', 'reverse.file_rev')
            copy_file('reverse.file_rev', one_block_ref_block_file)
            HD.close()
        HC.close()
    
        TT.write('\n')
    
        #对one_block_ref_block_file进行muscle比对
        #定义生成的muscle文件名称
        one_block_ref_block_file_muscle = ('one'+'-'+'block'+'-'+reference_strainID+'.'+str(subject_begin)+'-' +
                                           str(subject_end)+'.muscle')
        #运行muscle比对
        try:
            subprocess.run(['muscle', '-align', one_block_ref_block_file, '-output', one_block_ref_block_file_muscle],
                           check=True)
        #muscle命令操作异常处理
        except subprocess.CalledProcessError as e:
            print(f"muscle 命令执行失败: {e}")
        #print(one_block_ref_block_file_muscle)
    
        #**************由于部分block并不是所有基因组都具备，因此补充没被记录的基因组，并用'-----'表示序列长度***************************
        #**************如果有12个基因组，选择cov为7,那么会筛选出有些block比对上的参考基因组会没有序列，使用'-'代替**********************
        #下面就是记录所有比对上的基因组ID，以及muscle比对后的序列长度。再与所有基因组比较，发现哪些基因组没有比对上，就补充'-'***************
        #创建一个字典：ID为键，序列为值
        block_seq_id = {}
        #定义一个空ID名称
        block_one_seq_id = ''
        print('reading', one_block_ref_block_file_muscle, '...')
        #读取one_block_ref_block_file_muscle文件
        OR = open(one_block_ref_block_file_muscle, 'r')
        #按行读取文件
        muscle_content = OR.readlines()
        #遍历文件的每一行
        for z in range(0, len(muscle_content), 1):
            #如果'>'存在，记录基因组ID
            if '>' in muscle_content[z]:
                block_one_seq_id = muscle_content[z].split(' ')[0]
                block_one_seq_id = block_one_seq_id.split('>')[1]
                if '\n' in block_one_seq_id:
                    block_one_seq_id = block_one_seq_id.split('\n')[0]
    
            else:
                #如果字典中不存在这个键(基因组ID),记录该键的值为0
                if block_seq_id.get(block_one_seq_id, 'no_exist') == 'no_exist':
                    block_seq_id[block_one_seq_id] = 0
                #切掉该行的换行符(因为要记录序列的长度)
                muscle_content[z] = muscle_content[z].split('\n')[0]
                #记录序列长度：序列长度累加
                block_seq_id[block_one_seq_id] = block_seq_id[block_one_seq_id]+len(muscle_content[z])
        OR.close()
    
        #print(block_seq_id.keys())
        #print(genome_seq.keys())
        #遍历genome_seq中的键，该字典具有所有基因组的ID名称
        for id in genome_seq:
            #如果one_block_ref_block_file_muscle文件没有该基因组ID
            if block_seq_id.get(id, 'no_id') == 'no_id':
                #在one_block_ref_block_file_muscle文件后追加
                OG = open(one_block_ref_block_file_muscle, 'a+')
                #追加补充的'>基因组ID名称'以及换行符
                OG.write('>'+id+'\n')
                b_length = 0
                #遍历字典block_seq_id的键
                for b_id in block_seq_id:
                    #记录键值：序列长度
                    b_length = block_seq_id[b_id]
                    break
                #追加同样字符长度的'-'
                OG.write(b_length*'-'+'\n')
                OG.close()
    
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
