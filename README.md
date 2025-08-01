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
  * Install using a package manager.
    * Update the package repository.
      ```
      sudo apt update
      ```
    * Upgrade existing packages.
      ```
      sudo apt upgrade
      ```
    * Install Python3.
      ```
      sudo apt install python3
      ```
  * Install from source code (by compilation).
    * Download the Python source code: Download the source code package of the required version from the official Python website, for example, download Python 3.13.5.
      ```
      wget https://www.python.org/ftp/python/3.13.5/Python-3.13.5.tgz
      ```
    * Unzip and create an installation directory.
      ```
      tar -zxvf Python-3.13.5.tgz
      mkdir install_Python-3.13.5
      ```
    * Configure the compilation options (assuming the installation path is /home/install/python/install_Python-3.13.5/bin).
      ```
      cd Python-3.13.5
      ./configure --prefix=..../install_Python-3.13.5
      ```
    * Compile and install.
      ```
      make
      make install
      ```
    * Configure environment variables to enable global calls. Use the vi editor to modify the environment variable configuration file.
      ```
      vi ~/.bashrc
      ```
      Type "i" to enter "insert" mode (i.e., "editable mode"); add the configuration information at the end of the file.
      ```
      export PATH=/home/install/python/install_Python-3.13.5/bin:$PATH
      ```
      Press "Esc" to exit editing, then type ":wq" to save.
      Configure for immediate execution:
      ```
      source ~/.bashrc
      ```
* perl. Most Linux distributions come pre-installed with Perl.
  * Check if Perl is installed on the Linux system. If it is installed, the version information will be displayed.
    ```
    perl -v
    ```
  * Install using a package manager.
    ```
    sudo apt update
    sudo apt install perl
    ```
    
* Java. In Linux systems, Java typically exists in the form of OpenJDK or Oracle JDK. You can install OpenJDK via a package manager, for example:
  ```
  sudo apt update
  sudo apt update
  sudo apt install openjdk-11-jdk
  ```

* blast+ (available at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). Download ncbi-blast-2.17.0+-x64-linux.tar.gz.
  * Locate ncbi-blast-2.17.0+-x64-linux.tar.gz, right-click the mouse and select "Copy link", then use wget plus the copied link to download it:
    ```
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
    ```
  * Extract (the compressed file):
    ```
    tar -zxvf ncbi-blast-2.17.0+-x64-linux.tar.gz
    ```
  * Go to the directory and find the bin folder (this folder contains the running program, assuming the path to this folder is /home/install/blast/ncbi-blast-2.17.0+/bin), and configure the environment variables:
    * Edit the configuration file:
      ```
      vi ~/.bashrc
      ```
    * Add configuration:
      ```
      export PATH=/home/install/blast/ncbi-blast-2.17.0+/bin:$PATH
      ```
    * Configure for immediate execution:
      ```
      source ~/.bashrc
      ```

* Muscle v5.3 （available at https://www.drive5.com/muscle5/ ） <br>
  * Download the source code of Muscle v5.3 (this version is a binary file and does not need to be decompressed after download):
    ```
    wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3
    ```
  * Modify the file execution mode:
    ```
    chmod +x muscle-linux-x86.v5.3
    ```
  * Rename the file for easier calling:
    ```
    mv muscle-linux-x86.v5.3 muscle
    ```

* Gblocks 0.91b (available at https://www.biologiaevolutiva.org/jcastresana/Gblocks.html; https://slackbuilds.org/repository/14.2/academic/Gblocks/). <br>
  * Download Gblocks_Linux_0.91b.tar.Z and decompress it.
    ```
    wget https://www.biologiaevolutiva.org/jcastresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z
    tar -zxvf Gblocks_Linux64_0.91b.tar.Z
    ```
  * Go to the Gblocks_0.91b folder (for example, in /home/install/Gblocks/Gblocks_0.91b; you'll find the executable program Gblocks), and configure the environment variables to call Gblocks globally.
    * Use the vi editor to modify the environment variable configuration file.
      ```
      vi ~/.bashrc
      ```
    * Type "i" to enter "insert" mode (i.e., "editable mode"); add the configuration information at the end of the file.
      ```
      export PATH=/home/install/Gblocks/Gblocks_0.91b:$PATH
      ```
    * Press "Esc" to exit editing, then type ":wq" to save.
    * Configure for immediate execution:
      ```
      source ~/.bashrc
      ```
    * If you can call Gblocks in any directory and it can be recognized and used, it means the installation is complete.
* Jmodeltest v2.1.10 (available at https://github.com/ddarriba/jmodeltest2/releases).
  * Download Jmodeltest v2.1.10.
    ```
    wget https://github.com/ddarriba/jmodeltest2/archive/refs/tags/v2.1.10r20160303.tar.gz
    ```
  * Decompress jmodeltest-2.1.10.tar.gz.
    ```
    tar -zxvf jmodeltest-2.1.10.tar.gz
    ```
  * Configure the directory where jModelTest.jar is located.
    ```
    vi ~/.bashrc
    ```
  * Add configuration:
    ```
    export PATH=/home/install/jmodeltest/jmodeltest-2.1.10:$PATH
    ```
  * Configure for immediate execution:
    ```
    source ~/.bashrc
    ```
    
* raxml-ng v1.1.0 (available at https://github.com/amkozlov/raxml-ng/releases).
  * Download raxml-ng:
    ```
    wget https://github.com/amkozlov/raxml-ng/releases/download/1.1.0/raxml-ng_v1.1.0_linux_x86_64.zip
    ```
  * Decompress:
    ```
    unzip raxml-ng_v1.1.0_linux_x86_64.zip
    ```
  * Configure the environment (assuming the raxml-ng program is in the /home/install/raxml directory).
    ```
    vi ~/.bashrc
    ```
  * Add configuration:
    ```
    export PATH=/home/install/raxml:$PATH
    ```
  * Configure for immediate execution:
    ```
    source ~/.bashrc
    ```
* Forester.jar (available at https://github.com/cmzmasek/forester), Download it to the Windows system and double-click to run it.

How to use WMAF?
===============
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
=======
In this case, the mitogenomic data are from Saccharomyces species and downloaded from National Center for Biotechnology Information (NCBI: https://www.ncbi.nlm.nih.gov/). You can download it on this project (`00.integrate_13_Saccharomyces_genomes.txt`) and copy it under Linux system.<br>
*Check the file format.
  * Ensure that the file content is in Fasta format.
  * Files transferred from the Windows system to the Linux system need to have their format modified.
   ```
   less 00.integrate_13_Saccharomyces_genomes.txt |perl -ne 's/\r//g;print $_'>cc
   mv cc 00.integrate_13_Saccharomyces_genomes.txt
   ```
  * Ensure that the sequence file does not contain blank lines; if it does contain blank lines:
   ```
   less 00.integrate_13_Saccharomyces_genomes.txt |grep -v '^$' > cc
   mv cc 00.integrate_13_Saccharomyces_genomes.txt
   ```
* Run `whole-mitogenome_sequence_alignment.py` (using the S288c strain as the reference mitochondrial genome; set the minimum block length to 200; set the number of strains covering the reference genome to 13).
 ```
 python whole-mitogenome_sequence_alignment.py --all_genome 00.integrate_13_Saccharomyces_genomes.txt --reference_strainID S288c --block 200 --cov 13
 ```
* When the program finishes running, several block files will appear (for example: one-block-S288c.13774-13987.muscle). Use `integrate_aligned_blocks.py` to concatenate the alignment results of each block. The parameter `00.genomeID_list` required by the program has been automatically generated in the previous step.
  * Prepare the `block_path_list` file, which needs to be created manually. The command is as follows (the blocks need to be arranged in order):
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
  * Run the code below to concatenate each block. It is recommended to name the output file `08.integrate_blocks.txt`.
  ```
  python integrate_aligned_blocks.py --strainID 00.genomeID_list --block_path_list block_path_list --output 08.integrate_blocks.txt
  ```
* Use Gblocks to remove unreliable or highly variable sequence segments from the alignment results, thereby improving the accuracy of subsequent phylogenetic analyses (such as constructing phylogenetic trees). Name the generated file with the suffix "-gb" (for example: `08.integrate_blocks.txt-gb`).
* Modify the format of `08.integrate_blocks.txt-gb`.
  * The original format of the file is as follows:
    ```
    >CBS8840
    TATAAGTAAT AAATAAAGTT TATATAGTAT ATAAATATTA TATTATTATA TAATAAGGAA
    AATATTAGAT ATAAGAATAT GATGTTGGTT CAGATTAAGC GCTAAATAAG GACATGACAC
    ATGCGAATCA TACGTTTATT ATTGATAAGA TAATAAATAA GTGGTGTAAA CGTGAGTAAT
    ATATTAGGAA TAAAGAACTA TAGAATAAGC TAAATACTTA ATAAAAAATA TAAAAAATTA
    TATAAAAATA AAAGGAAAGT AATATTTATC TATAGTCAAG CCAATAATGG TTTAGGTAGT
    ```
  * Convert `08.integrate_blocks.txt-gb` to Fasta format:
    ```
    less 08.integrate_blocks.txt-gb |perl -ne 's/ //g;print "$_"'>cc
    mv cc 08.integrate_blocks.txt-gb
    ```
  * The modified file is as follows:
    ```
    >CBS8840
    TATAAGTAATAAATAAAGTTTATATAGTATATAAATATTATATTATTATATAATAAGGAA
    AATATTAGATATAAGAATATGATGTTGGTTCAGATTAAGCGCTAAATAAGGACATGACAC
    ATGCGAATCATACGTTTATTATTGATAAGATAATAAATAAGTGGTGTAAACGTGAGTAAT
    ATATTAGGAATAAAGAACTATAGAATAAGCTAAATACTTAATAAAAAATATAAAAAATTA
    TATAAAAATAAAAGGAAAGTAATATTTATCTATAGTCAAGCCAATAATGGTTTAGGTAGT
    ```
* Use Jmodeltest to predict the optimal model for constructing a phylogenetic tree (assuming the installation path of `JmodelTest.jar` is `/home/install/jmodeltest/jmodeltest-2.1.10/jModelTest.jar`):
  * Create a shell script file to run the program in the background using the nohup command.
  ```
  cat >09.run_jmodeltest.sh
  java -jar /home/install/jmodeltest/jmodeltest-2.1.10/jModelTest.jar -d 08.integrate_blocks.txt-gb -t ML -f -i -g 4 -AIC -BIC -a
  ^C
  ```
  * Run the program in the background.
  ```
  nohup sh 09.run_jmodeltest.sh > 09.run_jmodeltest.sh.log &
  ```
  * After the program finishes running, check the `09.run_jmodeltest.sh.log` file (the optimal model result is at the end of the file); the result shows that the optimal models are: AIC: GTR+I+G; BIC: GTR+I.
  ```
  ::Best Models::
  Model f(a) f(c) f(g) f(t) kappa titv Ra Rb Rc Rd Re Rf pInv gamma
  -----------------------------------------------------------------------------------------
  AIC GTR+I+G 0.37 0.10 0.13 0.39 0.00 0.00 1.514 5.527 8.292 1.413 10.585 1.000 0.44 0.17
  BIC GTR+I 0.37 0.10 0.13 0.39 0.00 0.00 1.480 5.384 7.994 1.414 10.348 1.000 0.83 N/A
  ```
* Construct a phylogenetic tree using the optimal model GTR+I+G. It is necessary to convert `08.integrate_blocks.txt-gb` to phylip format. Use `fasta2phylip.pl` for the conversion.
  ```
  fasta2phylip.pl 08.integrate_blocks.txt-gb > 08.integrate_blocks.txt-gb-phylip
  ```
* Use RAxML to construct a phylogenetic tree.
  ```
  cat >10.run_raxmal.sh
  raxml-ng --all --msa 08.integrate_blocks.txt-gb-phylip --threads 2 --model GTR+I+G4 --bs-trees 500
  ^C
  ```
  Run the program in the background.
  ```
  nohup sh 10.run_raxmal.sh > 10.run_raxmal.sh.log &
  ```
  Generate the file used for constructing the phylogenetic tree with the suffix ".support" (`08.integrate_blocks.txt-gb-phylip.raxml.support`).
* Use forester.jar to visualize the phylogenetic tree.
  * Click "`File`" in the navigation bar, then click "`Read Tree from File ...`".
  * Upload the "`.support`" file.
  * Click "`Options`", select "`Scale`" to display the phylogenetic tree scale; then click "`Select Color Scheme ...`" and choose "`Black & White`" to set the color background.
  * Click "`Tools`" and select "`Midpoint-Root`".
  * The phylogenetic tree will now be fully displayed.
  
* The phylogenetic tree after being refined with Inkscape is as follows:
![](https://github.com/cuihs-nl312/WMAF/blob/main/Example_data/08.inkescape.svg)
  
  

  

  

    
    


   
  




