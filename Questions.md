## 题目发布

<a name="content">目录</a>

- [项目说明](#introduction)
- [Q&A](#question-and-answer)
- [入门题](#for-beginer)
    - [1. 根据序列ID，提取目标序列](#for-beginer-1)
    - [2. 双端未匹配数据的重新匹配](#for-beginer-2)
- [进阶题](#for-user-with-middle-level)

- [挑战题](#for-veterans)

<a name="introduction"><h2>项目说明 [<sup>目录</sup>](#content)</h2></a>



<a name="question-and-answer"><h2>Q&A [<sup>目录</sup>](#content)</h2></a>

<a name="for-beginer"><h2>入门题 [<sup>目录</sup>](#content)</h2></a>

<a name="for-beginer-1"><h3>1. 根据序列ID，提取目标序列 [<sup>目录</sup>](#content)</h3></a>

根据序列ID，提取目标序列

给定：一个[fastq/fasta文件](./Attachments/R1.fastq)，以及另外一个保存目标序列的[id list文件](./Attachments/R1.interested.id)

任务：从fastq/fasta文件中，根据目标序列的ID，提取出目标序列保存为一个新的fastq/fasta文件

<a name="for-beginer-2"><h3>2. 双端未匹配数据的重新匹配 [<sup>目录</sup>](#content)</h3></a>

给定：一个样本的双端测序文件（Attachments文件夹下的R1.fastq和R2.fastq），且这双端Forward-end与Reverse-end同一行的序列并非如标准PE数据那样一一对应，即来源于同一个fragment

任务：只保留有双端序列的fragment，输出到处理后的双端Fastq文件，且让它在两个fastq文件的同一行一一对应

<a name="for-user-with-middle-level"><h2>进阶题 [<sup>目录</sup>](#content)</h2></a>



<a name="for-veterans"><h2>挑战题 [<sup>目录</sup>](#content)</h2></a>
