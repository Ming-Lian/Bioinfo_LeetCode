# -*- coding: utf-8 -*-
# @Time : 2019/8/7 10:26
# @Author : Zhongyi Hua
# @FileName: question3.py
# @Usage:LeeCode
# @Note:
# @E-mail: njbxhzy@hotmail.com

import os
import argparse
import csv

class MergeMatrix:
    def __init__(self, input_directory):
        self.dir = input_directory
        self.file = os.listdir(input_directory)
        self.pre_matrix = []
        self.matrix = []
        self.gene_set = set()
        self.sample_name = []

    def make_dict(self, filename):
        """
        Change every file to a dict, used for make_prematrix function to make pre_matrix
        :param filename: one sample expression file. First column should be geneID, second should be expression
        :return:{sample:name1,gene1:expression1..........}
        """
        with open(filename, "r") as f:
            lines = f.readlines()
            tmp_dict = {"sample": lines[0].strip().split("\t")[1]}
            for line in lines[1:]:
                tmp_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
            return tmp_dict

    def make_prematrix(self):
        """
        use every file dict to make expression list
        :return: [{sample:name1,gene1:expression1..........},
                  {sample:name2,gene1:expression1..........},
                  ......
                  {sample:name100,gene1:expression1..........}]
        """
        for _ in self.file:
            tmp_dict = self.make_dict(os.path.join(self.dir, _))
            self.pre_matrix.append(tmp_dict)
            self.sample_name.append(tmp_dict["sample"])
            tmp_set = set(tmp_dict.keys())
            tmp_set.remove("sample")
            self.gene_set |= tmp_set

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
        gene_index =dict()
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
    args = parser.parse_args()
    input_dir = args.input_directory
    output_path = args.output_path
    # main pipeline
    main_progress = MergeMatrix(input_dir)
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
