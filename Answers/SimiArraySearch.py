import sys
import argparse
import time

class Node :
    def __init__(self, initNumber, initArr):
        self.number = initNumber
        self.point = {}
        self.point[initArr] = 1
        self.next = None
    
    def getNum(self):
        return self.number
    
    def getPoint(self):
        return self.point
    
    def getNext(self):
        return self.next
        
    def setNum(self, newNum):
        self.number = newNum
    
    def addPoint(self, newArr):
        self.point[newArr] = self.point.get(newArr, 0) + 1
        
    def setNext(self, newNext):
        self.next = newNext

# 创建一个递增的有序链表类
class IncreaseList:
    def __init__(self):
        self.head = None
    
    def isEmpty(self):
        return self.head == None
        
    def printData(self):
        '''
        以Num1{ArrID1:n1;ArrID2:n2;...} Num2{...} ...的形式输出链表，用于检查链表是否创建正确
        '''
        current = self.head
        while current != None:
            print('%.5f{'%current.getNum(), end='')
            for ArrID in current.point.keys():
                print('%d:%d;'%(ArrID,current.point[ArrID]), end='')
            print('} ', end='')
            current = current.getNext()
        print('')
    
    def add(self, newNum, newArr):
        current = self.head
        previous = None
        found = False # 用于标记是否找到已有的值，只有没找到才需要添加新节点
        stop = False
        while current != None and not stop:
            # 找到已有的值，则在已有节点下添加一个ArrID->Count
            if current.getNum() == newNum:
                current.addPoint(newArr)
                found = True
                stop = True
            elif current.getNum() > newNum:
                stop = True
            else:
                previous = current
                current = current.getNext()
        if not found:
            temp = Node(newNum, newArr)
            if previous == None:
                temp.setNext(self.head)
                self.head = temp
            else:
                temp.setNext(current)
                previous.setNext(temp)

# 定义哈希链表类
class Hash2List:
    def __init__(self, resolution):
        '''
        初始化时，需要制定用于计算哈希值的分辨率
        '''
        self.hash = {}
        self.resolution = resolution
    
    def isEmpty(self):
        return len(self.hash.keys()) == 0
    
    def add(self, newNum, newArr):
        hashValue = float(str(newNum).split('.')[0] + '.' + str(newNum).split('.')[1][0:self.resolution])
        if not hashValue in self.hash:
            self.hash[hashValue] = IncreaseList()
        self.hash[hashValue].add(newNum, newArr)
    
    def searchSimilar(self, number, diff):
        start = number - diff
        end = number - diff
        start_hashValue = float(str(start).split('.')[0] + '.' + str(start).split('.')[1][0:self.resolution])
        end_hashValue = float(str(end).split('.')[0] + '.' + str(end).split('.')[1][0:self.resolution])
        found = False
        # 搜索起点和终点可能落在不同的区间
        for hashValue in {start_hashValue, end_hashValue}:
            if hashValue in self.hash:
                # 从该哈希值的链表的头开始查找
                current = self.hash[hashValue].head
                stop = False
                while True:
                    if current != None and not stop:
                        # 一旦超过搜索范围则停止
                        if current.getNum() - number > diff:
                            stop = True
                        else:
                            if abs(current.getNum() - number) <= diff:
                                found = True
                                for ArrID in current.point:
                                    yield (current.getNum(), ArrID)
                            current = current.getNext()
                    else:
                        break
        if not found:
            return

# 从文件中读取数据，并构造成倒排索引
def BuildHash2List_InvertIndex(filename, res=1):
    Hash2List_InvertIndex = Hash2List(res)
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                ArrID = int(line.strip('#'))
            else:
                Num = float(line.strip())
                Hash2List_InvertIndex.add(Num, ArrID)
    return Hash2List_InvertIndex

# 根据给定的倒排索引，对查询文件中的每个数组寻找索引文件的最相近数组
def CalcSimilar(Hash2List_InvertIndex, quryfile, outfile, diff):
    with open(outfile, 'w', encoding='utf-8') as F:
        with open(quryfile, 'r', encoding='utf-8') as f:
            SimiCount = {}
            Used = {}
            MostSimi = {}
            maxArr = None
            for line in f:
                if line.startswith('#'):
                    ArrID = int(line.strip('#'))
                    # 将上一个数组的计算结果写出，并且初始化该轮的相关变量
                    if len(SimiCount.keys()) > 0:
                        F.write('%d\t%d\n'%(maxArr, SimiCount[maxArr]))
                        SimiCount = {}
                        Used = {}
                        maxArr = None
                else:
                    Valid = False
                    MostSimi = {}
                    Num = float(line.strip())
                    currentSimi = Hash2List_InvertIndex.searchSimilar(Num, diff) # 这是一个迭代器
                    # 将所有落在[x-d, x+d]的所有值遍历一遍，找出每个数组的剩下的数值中最优取值
                    for currentNum, currentArrID in currentSimi:
                        # 先判断当前值是否已弹出，判断结果保存在变量Valid中
                        if currentArrID in Used:
                            if not currentNum in Used[currentArrID]:
                                Valid = True
                        else:
                            Valid = True
                        # 对当前最优数值进行更新
                        if Valid and abs(MostSimi.get(currentArrID, currentNum) - Num) >= abs(Num - currentNum):
                            MostSimi[currentArrID] = currentNum
                    # 进行打分和已使用数值弹出
                    for ArrID in MostSimi:
                        # 打分
                        SimiCount[ArrID] = SimiCount.get(ArrID, 0) + 1
                        if maxArr == None:
                            maxArr = ArrID
                        else:
                            if SimiCount[ArrID] > SimiCount[maxArr]:
                                maxArr = ArrID
                        # 数值弹出
                        if not ArrID in Used:
                            Used[ArrID] = set()
                        Used[ArrID].add(MostSimi[ArrID])
            F.write('%d\t%d\n'%(maxArr, SimiCount[maxArr]))

# 根据给定的Diff计算粗精度，只考虑Diff<1的情况
def getResolution(Diff):
    dot = False
    resolution = 0
    for char in str(Diff):
        if char == '.':
            dot = True
            continue
        if dot:
            if char != '0':
                break
            resolution += 1
    return resolution

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Buid for fast similar comparation between Target Arrays and Source Arrays')
    argparser.add_argument('-T', type=str, required=True, dest='Target', help='File for Target Arry')
    argparser.add_argument('-S', type=str, required=True, dest='Source', help='File for Source Arry')
    argparser.add_argument('-D', type=float, dest='Diff', default=0.01, help='Max diff for similar pairs [default=0.01]')
    argparser.add_argument('-O', type=str, dest='Out', default='R.txt', help='Output file')
    
    args = argparser.parse_args()
    
    start = time.time()
    # 根据Diff计算粗精度对应的小数点位数
    res = getResolution(args.Diff)
    Hash2List_InvertIndex = BuildHash2List_InvertIndex(args.Target, res)
    end_index = time.time()
    print('It takes %f s to build Invert Index for Target Arrys.'%(end_index - start))
    CalcSimilar(Hash2List_InvertIndex, args.Source, args.Out, args.Diff)
    end_compare = time.time()
    print('It takes %f s to compare each Source Arry against Target Arrys.'%(end_compare - end_index))
