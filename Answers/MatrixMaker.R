# 该脚本只写了一个函数，有两个参数：
# （1）子文件保存的文件夹路径
# （2）子文件的相同后缀，若后缀中有'.'，需要进行转义，写成'\\.'

MatrixMaker <- function(dir,suffix){
    files <- list.files(dir)
    features <- vector()    # 用于保存所有样本中出现的features
    samples <- vector()    # 用于保存所有样本的id
    # 获取所有样本中出现的unique克隆，与所有样本的id
    for (file in files) {
        if (grepl(suffix,file)){
            data <- read.table(paste(dir,file,sep="/"),header=F,sep="\t")
            features <- c(features,as.character(data$V1))
            samples <- c(samples,sub(suffix,"",file))
          }
    }
    # 创建用于保存所有样本整合的表达谱，以0填充
    features <- unique(features)
    dataMat <- matrix(rep(0,length(samples)*length(features)),nrow = length(features),ncol = length(samples))
    rownames(dataMat) <- features
    colnames(dataMat) <- samples
    # 逐一读入样本的count文件，对matrix中对应的元素进行赋值
    for (file in files) {
        if (grepl(suffix,file)){
            data <- read.table(paste(dir,file,sep="/"),header=F,sep="\t")
            sample <- as.character(sub(suffix,"",file))
            index <- match(data$V1,rownames(dataMat))
            dataMat[index,sample] <- data$V2
        }
    }
    dataMat
}
