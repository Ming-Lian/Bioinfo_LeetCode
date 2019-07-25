
![](./picture/BioLeetCode_Logo.png)

## 解题思路提示

<a name="content">目录</a>

- [入门题](#for-beginer)
    - [1. RNA序列翻译成蛋白质](#for-beginer-1)
    - [2. 获得反向互补序列](#for-beginer-2)
    - [3. 根据序列ID，提取目标序列](#for-beginer-3)
    - [4. 双端未匹配数据的重新匹配](#for-beginer-4)
    - [5. 将输入的大Fasta文件拆分成若干个小Fasta文件](#for-beginer-5)
- [进阶题](#for-user-with-middle-level)
    - [1. 从Fastq文件中随机抽样一定量的数据](#for-user-with-middle-level-1)
    - [2. 将输入的大矩阵文件按照列拆分成若干个sub-matrixs文件](#for-user-with-middle-level-2)
    - [3. 将若干个单样本的表达定量结果汇总成一个大矩阵，即expression profile matrix](#for-user-with-middle-level-3)

- [挑战题](#for-veterans)
    - [1. 分层Bootstrap抽样](#for-veterans-1)

<a name="for-beginer"><h2>入门题 [<sup>目录</sup>](#content)</h2></a>

<a name="for-beginer-1"><h3>1. RNA序列翻译成蛋白质 [<sup>目录</sup>](#content)</h3></a>

<a name="for-beginer-2"><h3>2. 获得反向互补序列 [<sup>目录</sup>](#content)</h3></a>

<a name="for-beginer-3"><h3>3. 根据序列ID，提取目标序列 [<sup>目录</sup>](#content)</h3></a>

（1）解析Fasta/Fastq文件

将给定的Fasta/Fastq文件逐行读入，解析ID行和序列行，若是Fastq文件，还需要解析质量行（它在Fastq文件的第4行）

注意：在Fastq文件中，每条记录的第3行为保留行，均为“+”字符，不包含任何信息，因此可以不需要进行解析，也不需要特意定义变量来保存这一行的信息

解析好之后，以哈希形式（在Perl中对键-值形式数据结构的叫法，在Python中成为字典）来存储解析后的数据，其中key为序列ID，value为序列字符串（若为Fastq文件，需要创建两个哈希，一个哈希保存：`ID -> Seq`，另一个哈希保存：`ID -> Qual`）

（2）根据ID list提取目标序列的记录

将给定的ID list文件逐行读入，若该ID所代表的序列在给定的Fasta/Fastq文件中存在，则在前面所创建的哈希变量中就能找到和它对应的键名(key)，则将其写出，否则给出提示：No such sequence!

示例代码，[点这里](./Answers/extractSeqFromFasta.pl)

<a name="for-beginer-4"><h3>4. 双端未匹配数据的重新匹配 [<sup>目录</sup>](#content)</h3></a>

(1) 解析Fastq文件

对双端的两个Fastq文件分别进行解析，保存成两组哈希，每组哈希有两个哈希变量，分别保存：`ID -> Seq` 和 `ID -> Qual`

解析方法与 [3. 根据序列ID，提取目标序列](#for-beginer-3)中的相同

(2) 找出双端能匹配的fragments将它们分别在输出的两个Fastq文件中在同一行写出

以两组哈希中的任意一组中的一个哈希的键进行遍历，例如以组一的`ID -> Seq`的键进行遍历，然后在另一组的`ID -> Seq`中查到对应的键是否存在，若存在，就说明双端是配对的，则进行输出

示例代码，[点这里](./Answers/PairsMate.pl)

<a name="for-beginer-5"><h3>5. 将输入的大Fasta文件拆分成若干个小Fasta文件 [<sup>目录</sup>](#content)</h3></a>

<a name="for-user-with-middle-level"><h2>进阶题 [<sup>目录</sup>](#content)</h2></a>

<a name="for-user-with-middle-level-1"><h3>1. 从Fastq文件中随机抽样一定量的数据 [<sup>目录</sup>](#content)</h3></a>
