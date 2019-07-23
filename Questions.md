## 题目发布

<a name="content">目录</a>

- [项目说明](#introduction)
- [Q&A](#question-and-answer)
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

<a name="introduction"><h2>项目说明 [<sup>目录</sup>](#content)</h2></a>

在IT圈中，LeetCode（中文名：力扣）是每一个程序猿修炼内功的闯关平台，几乎每一位程序员在顺利拿到IT大厂的offer之前，都被里面的各种数据结构与算法题虐过千百遍，但是正所谓不疯魔不成佛，每一个编程大牛都是在与一个一个编程任务和一个一个bug调试的斗争中逐渐成长起来的

与那些真正的动不动就写一个游戏，啪啪敲出成千上万行代码，而且还要考虑系统的稳健性和性能优越性的真程序猿相比，我们这些稍微懂点编程的生物狗，的确是一种非常辣鸡的存在

不过，也不要妄自菲薄，不是我们不行，而是我们工作生活中遇到的场景和他们不一样导致的，和他们相比，我们对系统的稳健性和性能优越性要求不高，基本上能把我要完成的数据分析任务搞定，那就是一个可行的编程解决方案，至于它的所谓时间复杂度与空间复杂度，不是需要我们考虑太多的东西，也没有所谓的用户或者产品经理给你提出比较尖酸刻薄的需求，没有压力，也就没有动力，所以慢慢的成了一条咸鱼

但是，只要你有一颗积极向上的心，只要还不想成为一个实验室的弱鸡搬砖人员，你的目标是星辰大海，既然外界不给咱们压力，那咱们就自己给自己造点压力——

<p align="center"><strong>朋友们，来接受编程的挑战吧</strong></p>

<a name="question-and-answer"><h2>Q&A [<sup>目录</sup>](#content)</h2></a>

<a name="for-beginer"><h2>入门题 [<sup>目录</sup>](#content)</h2></a>

<a name="for-beginer-1"><h3>1. RNA序列翻译成蛋白质 [<sup>目录</sup>](#content)</h3></a>

给定：氨基酸密码子表，和一个保存多条RNA序列的Fasta文件

任务：得到RNA序列翻译得到的氨基酸序列，保存为Fasta文件

<a name="for-beginer-2"><h3>2. 获得反向互补序列 [<sup>目录</sup>](#content)</h3></a>

给定：一条DNA序列

任务：得到这条DNA序列的反向互补序列

<a name="for-beginer-3"><h3>3. 根据序列ID，提取目标序列 [<sup>目录</sup>](#content)</h3></a>

根据序列ID，提取目标序列

给定：一个[fastq/fasta文件](./Attachments/R1.fastq)，以及另外一个保存目标序列的[id list文件](./Attachments/R1.interested.id)

任务：从fastq/fasta文件中，根据目标序列的ID，提取出目标序列保存为一个新的fastq/fasta文件

<a name="for-beginer-4"><h3>4. 双端未匹配数据的重新匹配 [<sup>目录</sup>](#content)</h3></a>

给定：一个样本的双端测序文件（Attachments文件夹下的R1.fastq和R2.fastq），且这双端Forward-end与Reverse-end同一行的序列并非如标准PE数据那样一一对应，即来源于同一个fragment

任务：只保留有双端序列的fragment，输出到处理后的双端Fastq文件，且让它在两个fastq文件的同一行一一对应

<a name="for-beginer-5"><h3>5. 将输入的大Fasta文件拆分成若干个小Fasta文件 [<sup>目录</sup>](#content)</h3></a>

给定：一个包含有大量序列的Fasta文件

任务：将该Fasta文件均匀拆分成若干个子Fasta文件，使得每个子文件中的序列条数相同（允许最后一个子文件的序列条数少一些），且提供两个拆分选项：

- 每个子文件的序列条数；
- 拆分得到的子文件的个数；

使得用户可以根据两个选项中的一个进行拆分

<a name="for-user-with-middle-level"><h2>进阶题 [<sup>目录</sup>](#content)</h2></a>

<a name="for-user-with-middle-level-1"><h3>1. 从Fastq文件中随机抽样一定量的数据 [<sup>目录</sup>](#content)</h3></a>


给定：双端FASTQ文件

任务：按照用户指定的数据量，即多少bp，从原始的双端FASTQ文件中随机抽取序列，要求双端FASTQ文件中的序列必须配对

<a name="for-user-with-middle-level-2"><h3>2. 将输入的大矩阵文件按照列拆分成若干个sub-matrixs文件 [<sup>目录</sup>](#content)</h3></a>

给定：给定一个大矩阵，行是features（基因或转录本），列是样本

任务：将该大矩阵按照列均匀拆分成若干个子sub-matrixes，使得每个子文件中的列数相同（允许最后一个子文件的列数少一些），且提供两个拆分选项：

- 每个子文件的列数；
- 拆分得到的子文件的个数；

使得用户可以根据两个选项中的一个进行拆分

<a name="for-user-with-middle-level-3"><h3>3. 将若干个单样本的表达定量结果汇总成一个大矩阵，即expression profile matrix [<sup>目录</sup>](#content)</h3></a>

给定：保存在同一个文件夹下的若干个单样本的表达定量文件，即有n行2列，每一行表示一个基因/转录本的定量，第一列是行名，即每每个基因/转录本的ID，第二列是样本

| GeneID | Sample_i |
|:---:|:---:|
| g1 | 12 |
| g2 | 44 |
| ... | ... |

任务：将将若干个单样本的表达定量结果汇总成一个大矩阵

要求：使用Perl、Python或R编写，从头实现，不允许调用现成的包或工具模块

<a name="for-veterans"><h2>挑战题 [<sup>目录</sup>](#content)</h2></a>

<a name="for-veterans-1"><h3>1. 分层Bootstrap抽样 [<sup>目录</sup>](#content)</h3></a>
