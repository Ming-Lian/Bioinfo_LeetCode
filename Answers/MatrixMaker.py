# -*- coding: utf-8 -*-
# @Time : 2019/8/7 10:26
# @Author : Zhongyi Hua
# @FileName: MatirxMaker.py
# @Usage:LeeCode
# @Note: Anotated By Lianm with Chinese
# @E-mail: njbxhzy@hotmail.com

import os
import argparse
import csv

# 将一个带有表头（分别为GeneId和SampleId）的两列的文件读入，保存为一个字典：
# {"sample" -> SampleId,
#  "gene1" -> expression1,
#  "gene2" -> expression2,
# ...}

def make_dict(filename):
    """
    Notice :This function fot file with head, to handle file without header, please use makedict2
    Change every file to a dict, used for make_prematrix function to make pre_matrix
    :param filename: one sample expression file. First column should be geneID, second should be expression
    :return:{sample:name1,gene1:expression1..........}
    """
    with open(filename, "r") as f:
        lines = f.readlines() # 一次性读入文件所有的行
        tmp_dict = {"sample": lines[0].strip().split("\t")[1]} # 由于文件的第一行的第一列保存的是样本Id，因此要从该行中将样本Id提取出来
        for line in lines[1:]:
            tmp_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
        return tmp_dict

# 该函数的实现目的与上面的函数相同，差别仅在于SampleId的提取来源不同，上面的函数是从文件的第一行中提取，而下面的函数是从文件名中提取的
def make_dict2(directory, filename):
    """
    Change every file to a dict, used for make_prematrix function to make pre_matrix
    :param filename: one sample expression file. First column should be geneID, second should be expression
            directory:self.dir
    :return:{sample:name1,gene1:expression1..........}
    """
    with open(os.path.join(directory, filename), "r") as f:
        lines = f.readlines()
        tmp_dict = {"sample": filename.strip().split(".")[0]} # 从文件名中提取SampleId
        for line in lines:
            tmp_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
        return tmp_dict


class MergeMatrix:
    def __init__(self, input_directory, header):
        self.dir = input_directory
        self.file = os.listdir(input_directory)
        self.pre_matrix = []
        self.matrix = []
        self.gene_set = set()
        self.sample_name = []
        self.header = header

    def make_prematrix(self):
        """
        use every file dict to make expression list
        :return: [{sample:name1,gene1:expression1..........},
                  {sample:name2,gene1:expression1..........},
                  ......
                  {sample:name100,gene1:expression1..........}]
        """
        if self.header is True:
            for _ in self.file:
                tmp_dict = make_dict(os.path.join(self.dir, _)) # 将当前文件解析成字典的形式：{"sample" -> SampleId, "gene1" -> expression1, "gene2" -> expression2, ...}
                self.pre_matrix.append(tmp_dict) # 将每个文件解析出来的字典变量进行追加保存
                self.sample_name.append(tmp_dict["sample"]) # 追加保存当前文件对应的SampleId，目的是得到总的SampleId List
                # 追加保存当前文件的GeneId List
                tmp_set = set(tmp_dict.keys())
                tmp_set.remove("sample")
                self.gene_set |= tmp_set
        else:
            for _ in self.file:
                tmp_dict = make_dict2(self.dir, _)
                self.pre_matrix.append(tmp_dict)
                self.sample_name.append(tmp_dict["sample"])
                tmp_set = set(tmp_dict.keys())
                tmp_set.remove("sample")
                self.gene_set |= tmp_set

    # 根据上面得到的GeneId List和SampleId List，构造初始值为0的矩阵，然后逐一修改填空
    def make_matrix(self):
        """
        get matrix from prematrix, the result in self.matrix
        :return: [[GeneID, sample1, sample2, sample3, sample4......],
                  [gene1,expression1,expression2.................. ],
                  ......
                  [gene100,expression1,expression2.................. ]]
        """
        self.matrix = [[0] * (self.sample_name.__len__()+1) for _ in range(self.gene_set.__len__()+1)]
        self.matrix[0] = self.sample_name
        self.matrix[0].insert(0, "GeneID")
        gene_list = list(self.gene_set)
        gene_list.insert(0, "zhanwei")
        sample_index = dict()
        gene_index = dict()
        for index, sample_name in enumerate(self.matrix[0]):
            sample_index[sample_name] = index
        for index, gene_name in enumerate(gene_list):
            gene_index[gene_name] = index
        for gene_name in gene_list[1:]:
            gene_row = [0] * (self.sample_name.__len__()+1)
            gene_row[0] = gene_name
            for gene_list in self.pre_matrix:
                sample_name = gene_list["sample"]
                gene_row[sample_index[sample_name]] = gene_list[gene_name]
            self.matrix[gene_index[gene_name]] = gene_row


if __name__ == '__main__':
    # parse cmd line
    parser = argparse.ArgumentParser(description="This script for merge gene matrix from mutiple files")
    parser.add_argument('-i', '--input_directory', required=True,
                        help='<filepath>  The filepath contain file')
    parser.add_argument('-o', '--output_path', nargs='?', const="stdout", type=str, help='<filepath>  output_path')
    parser.add_argument('--header', action="store_true", default=False, help='<bool>  if your data has header,please use --header',)
    args = parser.parse_args()
    input_dir = args.input_directory
    output_path = args.output_path
    headerstatus = args.header
    # main pipeline
    main_progress = MergeMatrix(input_dir, headerstatus)
    main_progress.make_prematrix()
    main_progress.make_matrix()
    if output_path == "stdout":
        for i in main_progress.matrix:
            print(i)
    else:
        csvFile = open(output_path, "w+")
        try:
            writer = csv.writer(csvFile)
            for i in range(main_progress.matrix.__len__()):
                writer.writerow(main_progress.matrix[i])
        finally:
            csvFile.close()
