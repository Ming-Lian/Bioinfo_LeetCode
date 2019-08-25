
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
    - [4. 利用Needleman–Wunsch 算法来编写一个简单的全局比对程序](#for-user-with-middle-level-4)
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

该任务允许用户指定两种拆分方式：指定子文件中的序列条数（为了方便下面的讨论说明，设为N）或者指定最终的输出子文件数（为了方便下面的讨论说明，设为M）——对于用户可能进行的不规范操作，即同时指定这两个条件，很可能这两个条件是相互矛盾的，此时只遵循前一个设定

不论用户最终实际上使用了哪种设定，都需要先解决一个问题：每一个子文件的序列条数n应该是多少条？

> - 如果用户已经直接指定了子文件的序列条数，那么就直接获取，不需要再去求解了，此时`n=N`；
>
> - 如果用户指定的是最终生成的子文件数，那么可以通过计算：`N=输入Fasta文件的中序列条数 / 子文件数= 输入Fasta文件的中序列条数 / M`，因为N可能不是整数，而序列条数是整数，肯定需要对N进行取整，那么如何取整呢？
>
>    输入Fasta文件的中序列条数，可以通过统计输入文件中以`'>'`起始的行数来确定
>
>    考虑一下上面可采用的商的取整的方案（向上取整`ceil` or 向下取整`floor`），以及不同的取整方案可能带来的影响：
>
>    - 若商不为整数，此时若向下取整，则意味着，实际上子文件中的序列条数比期望中的少一点，即n<N，则最终实际拆出的子文件会比用户指定的多至少一个，即`n<N，m>=M+1`；
>
>    - 若商不为整数，此时若向上取整，则意味着，实际上子文件中的序列条数比期望中的多一点，即n>N，那么最终实际上拆出的子文件数会小于等于用户指定的数量M，即m<=M
>
>        到底最后是等于还是小于，取决于之前的M-1个子文件中总共多写的序列条数是否大于n，即是否满足 `(N-n)(M-1) >= n`？当`(N-n)(M-1) < n`时，能保证正好拆出M个子文件，其中前M-1个子文件平均都有n条子序列，最后一个子文件的序列条数不满n条
>        
>        通过求解上面的不等式，可以算出n应该满足：`n > N(1- 1/M)`，而上面的向上取整已经限定了`n>N`，它满足`n > N(1- 1/M)`这个要求，**所以向上取整能够保证得到的子文件数正好等于M**

通过上面的推导论证，每一个子文件的序列条数n应该等于用户指定的N，或者`N=输入Fasta文件的中序列条数 / M`的向上取整

则剩下要做的就是逐行读入输入Fasta文件，通过对行首起始的`>`字符的识别来计数当前的Fasta序列的条数i，若i是能被n整除，则说明当前子文件刚好写满，此时要新起一个子文件

示例代码，[点这里](./Answers/splitFasta.pl)

<a name="for-user-with-middle-level"><h2>进阶题 [<sup>目录</sup>](#content)</h2></a>

<a name="for-user-with-middle-level-1"><h3>1. 从Fastq文件中随机抽样一定量的数据 [<sup>目录</sup>](#content)</h3></a>

<a name="for-user-with-middle-level-3"><h3>3. 将若干个单样本的表达定量结果汇总成一个大矩阵，即expression profile matrix [<sup>目录</sup>](#content)</h3></a>

本题有两个推荐的实现思路：

1. **通过构造双重哈希实现**

    （1）根据指定的文件夹（下面记作`dir`）和文件后缀（下面记作`pattern`），将文件路径为`dir/*pattern`的文件逐一读入，然后保存为双重哈希形式，即

    <p align="center">%H={feature → {sample → quant}}</p>

    同时记录下样本名的列表 @S，最后输出的矩阵的列所表示的样本的顺序，由样本名列表 @S 确定

    （2）然后，根据上面构造出来的双重哈希 %H 和样本名列表 @S，对每个 feature 逐行写出到输出文件中，若当前 feature 为 i，遍历样本名列表 @S，若当前样本名为 j，则 feature i 在样本 i 的取值记为 n<sub>ij</sub>，缺失值用0填充

    <p align="center">n<sub>ij</sub> = defined(H{i}{j}) ? H{i}{j} : 0</p>
    
    示例代码：[点这里（基于perl的实现）](./Answers/MatrixMaker.pl)

2. **先构造初始矩阵，然后再进行填充**

    （1）根据指定的文件夹（下面记作`dir`）和文件后缀（下面记作`pattern`），将文件路径为`dir/*pattern`的文件逐一读入，得到 unique feature list 和 sample list，为了方便后面的说明，分别设为变量 @feature_list 和 @sample_list

    （2）根据 @feature_list 和 @sample_list，构造初始矩阵，行数为 @feature_list 的长度，列数为 @sample_list 的长度，矩阵中的每一个元素的值都赋为0（以0为缺省值）

    矩阵的行与 @feature_list 对应，矩阵的列与 @sample_list 对应

    （3）再根据指定的文件夹（下面记作`dir`）和文件后缀（下面记作`pattern`），将文件路径为`dir/*pattern`的文件逐一读入：
    
    - 根据读入的文件的文件名，获取样本id，将它与 @sample_list 比较，从而确定对应的矩阵的列索引 j；
    
    - 根据读入文件中第一列的GeneId，将它与 @feature_list 比较，从而确定对应的矩阵的行索引 i；

    获得行索引 i 和列索引 j 之后，就修改矩阵中对应元素的值了
    
    注意：该方法适用于数据量比较小的情况，当数据量比较大时，推荐用第一种方法实现
    
    示例代码：[点这里（基于R的实现）](./Answers/MatrixMaker.R) [点这里（基于python的实现，该示例脚本由Zhongyi Hua同学提供）](./Answers/MatrixMaker.py) 

<a name="for-user-with-middle-level-4"><h3>4. 利用Needleman–Wunsch 算法来编写一个简单的全局比对程序[<sup>目录</sup>](#content)</h3></a>

（1）熟悉Needleman–Wunsch算法，可以在纸上利用二维矩阵画出最优的比对路线。

（2）然后你在纸上如何画的，你就如何写程序实现。如，先初始化一个二维矩阵，填写必要数据，按照公式，将矩阵中的每一小格都填好数据，其中每一小格需记录好，累积分数和最优来源。数据填完就开始回溯找出比对方案.

示例代码：[点这里](./Answers/max_similarity.py)
