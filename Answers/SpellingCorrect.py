import re, collections

# 把语料中的单词全部抽取出来，转成小写，并且去除单词中间的特殊符号
# 单词会成为字母序列，don't 就变成 don 和 t 了
def words(text): return re.findall('[a-z]+', text.lower()) 

# 训练一个概率模型，实际上就是数一数每个单词出现几次
def train(features):
    # 解决新词的平滑问题：从来没有过见过的新词一律假设出现过一次
    # 给任意的键设置一个默认值，使用一个匿名的 lambda:1 函数，设置默认值为 1
    model = collections.defaultdict(lambda: 1) 
    for f in features:
        model[f] += 1
    return model

NWORDS = train(words(file('big.txt').read())) # NWORDS[w] 存储了单词 w 在语料中出现了多少次


alphabet = 'abcdefghijklmnopqrstuvwxyz'

# 返回所有与单词 w 编辑距离为 1 的集合
# 个集合很大. 对于一个长度为 n 的单词，可能有n种删除，n-1中对换（临近位点的对换），26n 种 (实际上是 25n 种)替换，和 26(n+1) 种插入
def edits1(word):
    n = len(word)
    # 字符串索引是一个左闭右开的区间，即[a, b)
    return set([word[0:i]+word[i+1:] for i in range(n)] +                     # deletion
               [word[0:i]+word[i+1]+word[i]+word[i+2:] for i in range(n-1)] + # transposition
               [word[0:i]+c+word[i+1:] for i in range(n) for c in alphabet] + # alteration
               [word[0:i]+c+word[i:] for i in range(n+1) for c in alphabet])  # insertion

# 编辑距离为 2 的那些单词
# 一般讲拼写检查的文献宣称大约80-95%的拼写错误都是介于编译距离 1 以内。 然而下面我们看到，当我对于一个有270个拼写错误的语料做实验的时候，我发现只有76%的拼写错误是属于编辑距离为1的集合
# 这个事情很简单, 递归的来看, 就是把 edit1 函数再作用在 edit1 函数的返回集合的每一个元素上就行了. 因此, 我们定义函数 edit2:
# def edits2(word):
#    return set(e2 for e1 in edits1(word) for e2 in edits1(e1))
# 这个语句写起来很简单, 实际上背后是很庞大的计算量: 与 something 编辑距离为2的单词居然达到了 114,324 个.
# 不过编辑距离放宽到2以后, 我们基本上就能覆盖所有的情况了, 在270个样例中, 只有3个的编辑距离大于2

# 我们可以做一些小小的优化: 在这些编辑距离小于2的词中间, 只把那些正确的词作为候选词
# 只返回那些正确的并且与 w 编辑距离小于2 的词的集合
def known_edits2(word):
    # 现在, 在刚才的 something 例子中, known_edits2('something') 只能返回 3 个单词: 'smoothing', 'something' 和 'soothing', 而实际上所有编辑距离为 1 或者  2 的词一共有 114,324 个. 这个优化大约把速度提高了 10%
    return set(e2 for e1 in edits1(word) for e2 in edits1(e1) if e2 in NWORDS)

# 获取候选词库中在语料库中出现的单词
def known(words): return set(w for w in words if w in NWORDS)

# 建模P(w|c)。规则：编辑距离为1的正确单词比编辑距离为2的优先级高, 而编辑距离为0的正确单词优先级比编辑距离为1的高，即w0 > w1 > w2 > w'。这里不直接计算P(w|c)的绝对取值，只需要知道它们相对值大小
# 从一个候选集合中选取最大概率的. 实际上, 就是选取有最大 P(c) 值的那个. 所有的 P(c) 值都存储在 NWORDS 结构中
def correct(word):
    # 依据P(w|c)的构造规则，不直接计算P(w|c)的绝对取值，只需要知道它们相对值大小，放回优先级最高的那个级别的候选单词即可
    candidates = known([word]) or known(edits1(word)) or known_edits2(word) or [word]
    # 选取有最大 P(c) 值的那个候选单词
    return max(candidates, key=lambda w: NWORDS[w])
