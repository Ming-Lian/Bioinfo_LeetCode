## 解题思路提示

<a name="content">目录</a>

- [入门题](#for-beginer)
    - [1. RNA序列翻译成蛋白质](#for-beginer-1)
    - [2. 获得反向互补序列](#for-beginer-2)
    - [3. 根据序列ID，提取目标序列](#for-beginer-3)
    - [4. 双端未匹配数据的重新匹配](#for-beginer-4)
    - [5. 将输入的大Fasta文件拆分成若干个小Fasta文件](#for-beginer-5)
    - [6. 计算N50](#for-beginer-6)
    - [7. 计算测序深度(Coverage Depth)与覆盖度(Coverage Breadth)](#for-beginer-7)
    - [8. 生成长度为n的所有碱基序列](#for-beginer-8)
- [进阶题](#for-user-with-middle-level)
    - [1. 从Fastq文件中随机抽样一定量的数据](#for-user-with-middle-level-1)
    - [2. 将输入的大矩阵文件按照列拆分成若干个sub-matrixs文件](#for-user-with-middle-level-2)
    - [3. 将若干个单样本的表达定量结果汇总成一个大矩阵，即expression profile matrix](#for-user-with-middle-level-3)
- [挑战题](#for-veterans)
    - [1. 分层Bootstrap抽样](#for-veterans-1)
    - [2. 手写BWT](#for-veterans-2)
        - [2.1. Burrows-Wheeler Transformation](#for-veterans-2-1)
        - [2.2. BWT reverse transformation](#for-veterans-2-2)
    - [3. 手写BLAST](#for-veterans-3)
    - [4. 手写de Bruijn assembly](#for-veterans-4)
    - [5. 相似数组搜索](#for-veterans-5)
    - [6. 从头实现后缀树的序列比对：从树构建到序列比对](#for-veterans-6)

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

示例代码：[Perl版本](./Answers/extractSeqFromFasta.pl)

<a name="for-beginer-4"><h3>4. 双端未匹配数据的重新匹配 [<sup>目录</sup>](#content)</h3></a>

(1) 解析Fastq文件

对双端的两个Fastq文件分别进行解析，保存成两组哈希，每组哈希有两个哈希变量，分别保存：`ID -> Seq` 和 `ID -> Qual`

解析方法与 [3. 根据序列ID，提取目标序列](#for-beginer-3)中的相同

(2) 找出双端能匹配的fragments将它们分别在输出的两个Fastq文件中在同一行写出

以两组哈希中的任意一组中的一个哈希的键进行遍历，例如以组一的`ID -> Seq`的键进行遍历，然后在另一组的`ID -> Seq`中查到对应的键是否存在，若存在，就说明双端是配对的，则进行输出

示例代码：[Perl版本](./Answers/PairsMate.pl)

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

示例代码：[Perl版本](./Answers/splitFasta.pl)

<a name="for-beginer-6"><h3>6. 计算N50 [<sup>目录</sup>](#content)</h3></a>

可以分以下几步进行：

（1）读入fasta文件中的每条contig序列

（2）计算每条序列的长度，同时累加每条序列的长度得到总长度N

（3）按照序列长度从大到小遍历每条序列$n_i$，并逐一累加序列长度$L+=n_i$，一旦$L\ge 0.5N$，则结束遍历，则最后那条被访问的序列的长度即为N50

<p align='center'><img src=./picture/BioLeetCode_how2deal_easy-6.png /></p>

示例代码：[Perl版本](./Answers/calc_N50.pl)

<a name="for-beginer-7"><h3>7. 计算测序深度(Coverage Depth)与覆盖度(Coverage Breadth) [<sup>目录</sup>](#content)</h3></a>

<a name="for-beginer-8"><h3>8. 生成长度为n的所有碱基序列 [<sup>目录</sup>](#content)</h3></a>

可以按照下面的结构来构造序列组成：

<p align='center'><img src=./picture/BioLeetCode_how2deal_easy-8.png /></p>

对于长度的n的序列，总共有$4^n$种可能的序列，将其当作一个$4^n \times n$的矩阵来构造，逐列进行碱基组成的填充：

第一列，平均分成4个区块，，每个区块有$4^{n-1}$行，每个区块内的碱基组成相同，按照$A\to T \to C \to G$的顺序为每个区块选择碱基组成，对于其中的第$j$行，其碱基组成为由其所在的区块决定，其所在的块号为$b=\lceil i / 4^{n-1}\rceil$（$\lceil \cdot \rceil$表示向上取整），按照$A\to T \to C \to G$的顺序一个循环，其正好位于循环的第$b \% 4$个（$\%$表示取余运算）；

第二列，平均分成16个区块，，每个区块有$4^{n-2}$行，每个区块内的碱基组成相同，按照$A\to T \to C \to G$的顺序为每个区块选择碱基组成，对于其中的第$j$行，其碱基组成为由其所在的区块决定，其所在的块号为$b=\lceil i / 4^{n-2}\rceil$，按照$A\to T \to C \to G$的顺序一个循环，其正好位于循环的第$b \% 4$个；

...

第$i$列，平均分成$4^i$个区块，，每个区块有$4^{n-i}$行，每个区块内的碱基组成相同，按照$A\to T \to C \to G$的顺序为每个区块选择碱基组成，对于其中的第$j$行，其碱基组成为由其所在的区块决定，其所在的块号为$b=\lceil i / 4^{n-j}\rceil$，按照$A\to T \to C \to G$的顺序一个循环，其正好位于循环的第$b \% 4$个

示例代码：[Perl版本](./Answers/basePermutation.pl)

<a name="for-user-with-middle-level"><h2>进阶题 [<sup>目录</sup>](#content)</h2></a>

<a name="for-user-with-middle-level-1"><h3>1. 从Fastq文件中随机抽样一定量的数据 [<sup>目录</sup>](#content)</h3></a>

首先，需要对要完成的任务作一个简单的数学推导：

假设原始给定的输入序列条数为N（对于双端配对的测序数据，即为单端数据的序列条数），目标抽样数据量为G(bp)，已知测序长度为l x 2（假设序列等长），则可知目标抽样序列条数为：

$$n=\frac{G}{2l}$$

则，现在问题变成了要从N条序列中随机抽出n条序列，可以把这个问题等效于进行N次伯努利实验，随机变量为“X=得到‘是’的次数”，每次得到“是”的概率为$\frac{n}{N}$，则此时随机变量服从二项分布，即$X\sim b(N,\frac{n}{N})$，且若实验结果为“是”，则抽中这条序列，那么最终被抽中的序列条数为二项分布的期望：

$$E(X)=N . \frac{n}{N}=n$$

那么，如何编程模拟这个随机实验过程呢？

随机生成0~1之间的随机数$x$（均匀分布），将它与$\frac{n}{N}$比较，若$x\le\frac{n}{N}$则实验结果为“是”，则结果为“是”的概率正好为$\frac{n}{N}$

示例代码：[Perl版本](./Answers/randExtract.pl)

<a name="for-user-with-middle-level-3"><h3>3. 将若干个单样本的表达定量结果汇总成一个大矩阵，即expression profile matrix [<sup>目录</sup>](#content)</h3></a>

本题有两个推荐的实现思路：

1. **通过构造双重哈希实现**

    （1）根据指定的文件夹（下面记作`dir`）和文件后缀（下面记作`pattern`），将文件路径为`dir/*pattern`的文件逐一读入，然后保存为双重哈希形式，即

    <p align="center">%H={feature → {sample → quant}}</p>

    同时记录下样本名的列表 @S，最后输出的矩阵的列所表示的样本的顺序，由样本名列表 @S 确定

    （2）然后，根据上面构造出来的双重哈希 %H 和样本名列表 @S，对每个 feature 逐行写出到输出文件中，若当前 feature 为 i，遍历样本名列表 @S，若当前样本名为 j，则 feature i 在样本 i 的取值记为 n<sub>ij</sub>，缺失值用0填充

    <p align="center">n<sub>ij</sub> = defined(H{i}{j}) ? H{i}{j} : 0</p>

    示例代码：[Perl版本](./Answers/MatrixMaker.pl)

2. **先构造初始矩阵，然后再进行填充**

    （1）根据指定的文件夹（下面记作`dir`）和文件后缀（下面记作`pattern`），将文件路径为`dir/*pattern`的文件逐一读入，得到 unique feature list 和 sample list，为了方便后面的说明，分别设为变量 @feature_list 和 @sample_list

    （2）根据 @feature_list 和 @sample_list，构造初始矩阵，行数为 @feature_list 的长度，列数为 @sample_list 的长度，矩阵中的每一个元素的值都赋为0（以0为缺省值）

    矩阵的行与 @feature_list 对应，矩阵的列与 @sample_list 对应

    （3）再根据指定的文件夹（下面记作`dir`）和文件后缀（下面记作`pattern`），将文件路径为`dir/*pattern`的文件逐一读入：
    
    - 根据读入的文件的文件名，获取样本id，将它与 @sample_list 比较，从而确定对应的矩阵的列索引 j；
    
    - 根据读入文件中第一列的GeneId，将它与 @feature_list 比较，从而确定对应的矩阵的行索引 i；

    获得行索引 i 和列索引 j 之后，就修改矩阵中对应元素的值了

    示例代码: [R版本](./Answers/MatrixMaker.R) [Python版本](./Answers/MatrixMaker.py)

<a name="for-veterans"><h2>挑战题 [<sup>目录</sup>](#content)</h2></a>

<a name="for-veterans-2"><h3>2. 手写BWT [<sup>目录</sup>](#content)</h3></a>

<a name="for-veterans-2-1"><h4>2.1. Burrows-Wheeler Transformation [<sup>目录</sup>](#content)</h4></a>

其实，这一小题里的任务，以及算法的基本操作过程已经很清楚了，关键是采用什么样的数据结构去实现它

最简单的选择就是数组（以下以Perl语言的逻辑进行说明）：

> 例如：给定的序列为 `TCATC`
> 
> 首先，将其打散成一个字符一个元素的数组，变成`@seq=('T','C','A','T','C')`，并在数组最后追加上一个终止标识符`$`，从而变成`@seq=('T','C','A','T','C','\$')`（注意这个终止标识符属于特殊字符，需要进行转义）
> 
> 随后，将这个数组当作一个队列使用，对队尾元素先出队，然后再将其在队首入队，例如：
> 
> ```perl
> $tail = pop @seq;
> unshift @seq, $tail;
> ```
> 
> 重复执行这样的操作直到从队尾取出一个元素后，新的队尾是终止标识符时截至
> 
> 每执行一次，就将新队形的数组合并成字符串，追加保存到一个新数组中，例如`@BWT`：
> 
> ```perl
>  $queue = join "",@seq;
>  push @BWT, $queue;
> ```
> 最后，将`@BWT`中的字符串元素按照字母顺序逐一进行处理：提取最后一个字母

示例代码：[Perl版本](./Answers/BWT_index.pl)

<a name="for-veterans-2-2"><h4>2.2. BWT reverse transformation [<sup>目录</sup>](#content)</h4></a>

(1)方法一

在任务说明中，已经列出了该方法实现BWT逆变换的过程：

> 初始条件为：当前碱基位置index=0，即为FC的第一行，且当前碱基组成为base=\$，当前参考序列组成为T=\$
> 
> 循环执行下面的操作，直至当前碱基组成再一次为\$：
> 
> （1）回溯：获取LC同处于index行的碱基组成c，更新当前参考序列T=cT
> 
> （2）LC到FC的定位：获取LC处于index行的碱基c在FC对应的碱基的位置i，更新当前碱基位置index=i

<p align='center'><img src=./picture/BioLeetCode_issue_Hard_2-2-1.png height=300/></p>

那么现在的问题是：该如何分别实现上面的**回溯**和**LC到FC的定位**呢？

对于**回溯**问题，很好解决，我只需要用到基于BWT得到的Last Column（下面为了方便讨论，将其称为BWT数组），则根据已知条件里给定的当前碱基位置index，则可以直接回溯得到LC同处于index行的碱基组成c=BWT[index]

对于**LC到FC的定位**问题，由于FC列中的碱基组成是严格按照碱基优先级递增排列的，因此，只需要知道在参考序列T中碱基优先级小于当前碱基c的所有可能的碱基类型的计数总和Pre(c)，再加上当前碱基是同类碱基的计数次数OCC(c)，如下图：

<p align='center'><img src=./picture/BioLeetCode_how2deal_Hard_2-2-1.png height=300/></p>

则LC中索引位置为index的碱基c，它再FC中对应的同一碱基的索引位置为：

$$\mathrm{LF(c,index)=Pre(c) + OCC(c)}$$

因此，为了能够实现以上操作，我们需要额外增加以下两个索引信息：

> - **Pre散列**：键为各种可能的碱基类型，键值为参考序列中碱基优先级小于对应碱基的碱基总数
> - **OCC数组**：LC对应行碱基，在LC该行及该行之前出现的次数（0-base）

<p align='center'><img src=./picture/BioLeetCode_how2deal_Hard_2-2-2.png height=400/></p>

其中BWT数组和OCC数组是FM-index的组成部分，除此之外还包含SA数组（在BWT逆变换中还暂时用不到它，所以在该小节未设置该变量，它在下一小节BWT search中会用到），详细的FM-index，见下方说明：

> 后缀数组（suffix array）的FM-index：
>
> <p align='center'><img src=./picture/BioLeetCode_how2deal_Hard_2-2-3.png height=400/></p>
>
> 以上图的BWT matrix为例，后缀数组由下面三个数组组成：
>
> - **SA数组**：LC对应行碱基在原始序列T中的位置（1-base）
> - **OCC数组**：LC对应行碱基，在LC该行及该行之前出现的次数（0-base），例如对于LC的第5行（0-base）的碱基为a，在LC的第5行及第5行之前，a碱基只出现2次，故OCC[5]=1
> - **BWT数组**：即BWT output，上图的给的例子中即为`gc$aaac`

有了上面定义的三个参考序列的索引信息（BWT数组、OCC数组和Pre散列）后，BWT逆变换算法过程，就可以表示为以下形式：

> 给定：BWT数组、OCC数组和Pre散列
> 
> 初始条件为：当前碱基位置index=0，即为FC的第一行，且当前碱基组成为base=\$，当前参考序列组成为T=\$
> 
> 循环执行下面的操作，直至当前碱基组成再一次为\$：
> 
> （1）回溯：获取LC同处于index行的碱基组成c，即c=BWT[index]，更新当前参考序列T=cT
> 
> （2）LC到FC的定位：获取LC处于index行的碱基c在FC对应的碱基的位置i，即i=Pre[c] + OCC[index]，更新当前碱基位置index=i

示例代码：[Perl版本](./Answers/BWT_revTrans_V1.pl)

(2)方法二

<p align='center'><img src=./picture/BioLeetCode_issue_Hard_2-2-2.png height=300/></p>

可以通过设置两个数组来实现：

> - BWT数组：保存BWT转换后的BWT索引，即@BWT=('C', 'C', 'T', 'T', 'A', '$')
> - Array数组：保存已经得到的部分BWT矩阵，每一个元素为部分BWT矩阵的一行

按照以下流程来完成完整的BWT矩阵的构建：

> - 初始状态：
>
>    @BWT=('C', 'C', 'T', 'T', 'A', '$')
>
>    @Array为空
>
> - 循环执行以下操作，直到Array的列数达到@BWT的长度：
> 
>   （1）@Array**之前**追加一列，为@BWT；
> 
>   （2）对@Array的行按照字母表顺序进行排序，得到重排后的BWT矩阵，以此来更新@Array；

示例代码：[Perl版本](./Answers/BWT_revTrans_V2.pl)

<a name="for-veterans-5"><h3>5. 相似数组搜索 [<sup>目录</sup>](#content)</h3></a>

算法实现的逻辑为：

> 对应一个给定的查询数组（来自于`S.txt`），逐一对每个数值进行比较；
>
> 对应当前目标数值x，快速在`P.txt`文件中，快速找到落在 `[x-d, x+d ]` (其中d表示定义的最远相似距离) 的数值集合 `N={n1, n2, …}`，以及集合N中每个数值对应的`P.txt`中数值来源`A={a1, a2, …}`， 然后对找到的数值计分，计一分；
>
> 当一个数组遍历比较结束后，找出得分最高的`P.txt`中的数组

那么这个问题的关键就在于：**怎么快速找到落在`[x-d, x+d ]` 的数值集合`N={n1, n2, …}`，以及集合N中每个数值对应的`P.txt`中数值来源`A={a1, a2, …}` ？**

关键在于**对 `P.txt` 的高效索引**，主要采用了**倒排索引 + 哈希链表**的方法

> 1. **采用了数字到数组的倒排索引**
>
>       对 `P.txt` 文件中的每个数值，建立数值指向其来源数组的倒排索引，一旦定位到了这个数值就可以马上得到这个数值来源于哪个数组
>
> 2.	**以哈希链表形式来组织`P.txt`中的数值**（即倒排索引的字典）
>
>       目标搜索的数值范围是 `[x-d, x+d]`，但是一开始的搜索起点 `start` 肯定要 `start <=x-d` ，最后的搜索终点end肯定要保证 `end >= x+d`，即实际的搜索范围 `[start, end]` 肯定要包含 `[x-d, x+d]`，而为了提高搜索效率，要尽量缩小搜索范围，减少不必要的搜索，对此我采用哈希链表的方法：
>
>           以d的从左到右第一位非0有效数字的更高一位作为一个粗精度，例如d=0.01，则我们以0.1作为一个粗精度，将落在[X.1, X.2 )
>           范围的数值存在一个递增链表里，并用X.1这个值作为哈希的键，指向这个递增链表，则根据键，就可以知道其所存储的链表的数值区间
>
>       	对于一个给定数值x，它的目标搜索区间是[x-d, x+d]，这个区间要么落在一个咱们上面定义的区间内，要么横跨相邻的两个区间，
>           则在一个或两个区间内搜索就可以了

以上是算法性能提升的主要实现原理，除此之外，还在以下几个小细节处，进行了优化：

- 得到最高得分的数组

    当对一个数组中的每个数值进行相似比较计分后，最后需要从计分结果里把得分最高的那个数组找出来，一种最简单粗暴的方法是：最后遍历整个得分结果，找出最高分的那个数组

    但是这样做太低效，我在这里设置了一个变量来记录进行每一个数值的比较打分后，当前得分最高的数组，若有更高分，则对它进行更新

    这样一个查询数组里的每个数值比较打分结束后，直接就得到了得分最高的数组

示例代码：[Python版本](./Answers/SimiArraySearch.py)

补充学习资料：

> - 倒排索引：[什么是倒排索引？](https://blog.csdn.net/starzhou/article/details/87519973)
>
> - 链表：[Python实现链表](https://www.cnblogs.com/wangxiayun/p/8358991.html)
