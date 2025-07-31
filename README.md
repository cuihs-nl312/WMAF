WMAF User Guide
===============

Overview
=============
WMAF is a pipeline for whole-mitogenome sequence alignment in fungi, implemented in Python. It utilizes pairwise sequence alignment to retrieve conserved regions (shared blocks) of all mitogenomes, which are then used to construct the phylogenetic tree. 

Keys to successful run WMAF pipeline
====================================

Ubuntu running
--------------
Windows users can use the Windows Subsystem for Linux (WSL) to run WMAF. And you can download Ubuntu to run WMAF.

Requirement 
-----------
* python3 or higher is required <br>
  * 使用包管理器安装：
    * 更新软件包存储库，运行命令：
      ```
      sudo apt update
      ```
    * 升级现有软件包，运行命令：
      ```
      sudo apt upgrade
      ```
    * 安装 Python 3，运行命令：
      ```
      sudo apt install python3
      ```
  * 源码编译安装：
    * 下载 Python 源码：从Python 官方网站下载所需版本的源码包，例如下载 Python 3.13.5
      ```
      wget https://www.python.org/ftp/python/3.13.5/Python-3.13.5.tgz
      ```
    * 解压并创建安装目录：
      ```
      tar -zxvf Python-3.13.5.tgz
      mkdir install_Python-3.13.5
      ```
    * 配置编译选项，（假设安装路径是/home/install/python/install_Python-3.13.5/bin）：
      ```
      cd Python-3.13.5
      ./configure --prefix=..../install_Python-3.13.5
      ```
    * 编译并安装：
      ```
      make
      make install
      ```
    * 配置环境变量即可全局调用：
      ```
      vi ~/.bashrc
      ```
      添加配置：
      ```
      export PATH=/home/install/python/install_Python-3.13.5/bin:$PATH
      ```
      配置立即执行：
      ```
      source ~/.bashrc
      ```
 
* Muscle v5.3 （available at https://www.drive5.com/muscle5/ ） <br>
  * 下载Muscle v5.3源码(此版本为二进制文件，下载后不需要解压缩)：
    ```
    wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3
    ```
  * 修改文件执行方式：
    ```
    chmod +x muscle-linux-x86.v5.3
    ```
  * 修改文件名称便于调用：
    ```
    mv muscle-linux-x86.v5.3 muscle
    ```

How to use WMAF?
----------------
You need to download mitogenomic data and write into a fasta format file (Note: Please rename the file so that it ends with .txt).WMAF has been integrated into whole-mitogenome_sequence_alignment.py for sequence alignment at the whole-genome level. You can run whole-mitogenome_sequence_alignment.py using the command below.<br>
```
python whole-mitogenome_sequence_alignment.py --all_genome ALL_GENOME --reference_strainID REFERENCE_STRAINID --block BLOCK --cov COV
```
If you need assistance with parameter usage, you can use the -h option to obtain detailed information.
```
python whole-mitogenome_sequence_alignment.py -h
```
Detailed parameter information is provided in the table below.

Parameter  | detailed information|
--------- | --------|
all_genome  | a fasta file containing all genomes; Please modify the file name to end with .txt. |
reference_strainID  | reference strain ID, such as 'DRR332962' |
block  | The size of the smallest block, such as '200' |
cov  | The smallest number of genomes aligned to the reference genome |

Example
-------
In this case, the mitogenomic data are from Saccharomyces species and downloaded from National Center for Biotechnology Information (NCBI: https://www.ncbi.nlm.nih.gov/). You can download it on this project (`00.integrate_13_Saccharomyces_genomes.txt`) and copy it under Linux system.<br>
* 数据格式检查：
  * 确保文件内容为Fasta格式。
  * windows系统转换到Linux系统需要修改格式：
   ```
   less 00.integrate_13_Saccharomyces_genomes.txt |perl -ne 's/\r//g;print $_'>cc
   mv cc 00.integrate_13_Saccharomyces_genomes.txt
   ```
  * 确保序列文件不包含空行，如果包含空行：
   ```
   less 00.integrate_13_Saccharomyces_genomes.txt |grep -v '^$' > cc
   mv cc 00.integrate_13_Saccharomyces_genomes.txt
   ```
* 运行`whole-mitogenome_sequence_alignment.py` (以S288c菌株为参考线粒体基因组；最小block长度设为200；比对上参考基因组的菌株覆盖数为13)。
 ```
 python whole-mitogenome_sequence_alignment.py --all_genome 00.integrate_13_Saccharomyces_genomes.txt --reference_strainID S288c --block 200 --cov 13
 ```
* 程序运行结束，会出现若干block文件(例如：one-block-S288c.13774-13987.muscle）。使用`integrate_aligned_blocks.py`串联各blocks的比对结果。程序所需参数`00.genomeID_list`已经由上一步自动生成。
  * 准备`block_path_list`文件，该文件需要手动创建，命令如下（需要将block按顺序排列）：
  ```
  cat >block_path_list
  one-block-S288c.6442-6986.muscle
  one-block-S288c.6987-7986.muscle
  one-block-S288c.7987-8196.muscle
  one-block-S288c.13774-13987.muscle
  one-block-S288c.20652-20974.muscle
  one-block-S288c.26231-26707.muscle
  one-block-S288c.28433-29273.muscle
  one-block-S288c.36508-36935.muscle
  one-block-S288c.40842-41093.muscle
  one-block-S288c.43297-43649.muscle
  one-block-S288c.46666-46988.muscle
  one-block-S288c.58141-59140.muscle
  one-block-S288c.59141-60140.muscle
  one-block-S288c.60141-60396.muscle
  one-block-S288c.61867-62355.muscle
  one-block-S288c.66086-66309.muscle
  one-block-S288c.70091-70320.muscle
  one-block-S288c.73697-74534.muscle
  one-block-S288c.77342-77561.muscle
  one-block-S288c.79164-80023.muscle
  ^C
  ```
  * 运行下方代码，串联各block，输出文件建议命名为08.integrate_blocks.txt。
  ```
  python integrate_aligned_blocks.py --strainID 00.genomeID_list --block_path_list block_path_list --output 08.integrate_blocks.txt
  ```
* 使用Gblocks去除比对结果中不可靠或高变区的序列片段，从而提高后续系统发育分析（如构建进化树）的准确性。生成文件以-gb表示（如：08.integrate_blocks.txt-gb）。
* 修改08.integrate_blocks.txt-gb格式。
  * 文件原格式，如：
    ```
    >CBS8840
    TATAAGTAAT AAATAAAGTT TATATAGTAT ATAAATATTA TATTATTATA TAATAAGGAA
    AATATTAGAT ATAAGAATAT GATGTTGGTT CAGATTAAGC GCTAAATAAG GACATGACAC
    ATGCGAATCA TACGTTTATT ATTGATAAGA TAATAAATAA GTGGTGTAAA CGTGAGTAAT
    ATATTAGGAA TAAAGAACTA TAGAATAAGC TAAATACTTA ATAAAAAATA TAAAAAATTA
    TATAAAAATA AAAGGAAAGT AATATTTATC TATAGTCAAG CCAATAATGG TTTAGGTAGT
    ```
  * 将格式修改成Fasta格式：
    ```
    less 08.integrate_blocks.txt-gb |perl -ne 's/ //g;print "$_"'>cc
    mv cc 08.integrate_blocks.txt-gb
    ```
  * 修改后的文件格式如下：
    ```
    >CBS8840
    TATAAGTAATAAATAAAGTTTATATAGTATATAAATATTATATTATTATATAATAAGGAA
    AATATTAGATATAAGAATATGATGTTGGTTCAGATTAAGCGCTAAATAAGGACATGACAC
    ATGCGAATCATACGTTTATTATTGATAAGATAATAAATAAGTGGTGTAAACGTGAGTAAT
    ATATTAGGAATAAAGAACTATAGAATAAGCTAAATACTTAATAAAAAATATAAAAAATTA
    TATAAAAATAAAAGGAAAGTAATATTTATCTATAGTCAAGCCAATAATGGTTTAGGTAGT
    ```
* 使用Jmodeltest预测构建系统发育树的最优模型(假设安装Jmodeltest.jar安装路径为/home/install/jmodeltest/jmodeltest-2.1.10/jModelTest.jar)：
  * 创建shell文件：
  ```
  cat >09.run_jmodeltest.sh
  java -jar /home/install/jmodeltest/jmodeltest-2.1.10/jModelTest.jar -d 08.integrate_blocks.txt-gb -t ML -f -i -g 4 -AIC -BIC -a
  ^C
  ```
  * 后台运行程序：
  ```
  nohup sh 09.run_jmodeltest.sh > 09.run_jmodeltest.sh.log &
  ```
  * 程序运行结束后，查看`09.run_jmodeltest.sh.log`文件（最优模型结果在文件末尾）;结果显示最优模型为AIC：GTR+I+G；BIC：GTR+I。
  ```
  ::Best Models::
  Model f(a) f(c) f(g) f(t) kappa titv Ra Rb Rc Rd Re Rf pInv gamma
  -----------------------------------------------------------------------------------------
  AIC GTR+I+G 0.37 0.10 0.13 0.39 0.00 0.00 1.514 5.527 8.292 1.413 10.585 1.000 0.44 0.17
  BIC GTR+I 0.37 0.10 0.13 0.39 0.00 0.00 1.480 5.384 7.994 1.414 10.348 1.000 0.83 N/A
  ```
* 以最优模型GTR+I+G，构建系统发育树。需要将`08.integrate_blocks.txt-gb`转换成phylip格式。使用`fasta2phylip.pl`进行修改。
  ```
  fasta2phylip.pl 08.integrate_blocks.txt-gb > 08.integrate_blocks.txt-gb-phylip
  ```
* 使用Raxmal构建系统发育树。
  ```
  cat >10.run_raxmal.sh
  raxml-ng --all --msa 08.integrate_blocks.txt-gb-phylip --threads 2 --model GTR+I+G4 --bs-trees 500
  ^C
  ```
  后台运行程序：
  ```
  nohup sh 10.run_raxmal.sh > 10.run_raxmal.sh.log &
  ```
  生成的构建发育树的文件以.support结尾（`08.integrate_blocks.txt-gb-phylip.raxml.support`）。
* 可视化系统发育树（使用forester.jar），经Inkscape美化后的发育树如下图：
![](https://github.com/cuihs-nl312/WMAF/blob/main/Example_data/08.inkescape.svg)
  
  

  

  

    
    


   
  




