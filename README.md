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
You need to download mitogenomic data and write into a fasta format file (Note: Please rename the file so that it ends with .txt.)
WMAF has been integrated into whole-mitogenome_sequence_alignment.py for sequence alignment at the whole-genome level. You can run whole-mitogenome_sequence_alignment.py using the command below.<br>
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
In this case, the mitogenomic data are from Saccharomyces species and downloaded from National Center for Biotechnology Information (NCBI: https://www.ncbi.nlm.nih.gov/). You can download it on this project (`example_data`).



