library(xbreed)

setwd('/Users/zhezhang/Desktop/icloud/0-tongbu/6-teacher/1-book/CCPS2')

## 第一步：创建历史群体

# 基因组（3条染色体）
genome<-data.frame(matrix(NA, nrow=3, ncol=6))
names(genome)<-c("chr","len","nmrk","mpos","nqtl","qpos")
genome$chr<-c(1:3)
genome$len<-c(80,60,50)
genome$nmrk<-c(500,1000,250)
genome$mpos<-c('rnd','rnd','rnd')
genome$nqtl<-c(40,50,45)
genome$qpos<-c('rnd','rnd','rnd')
genome

# 500个世代，最终形成200个个体组成的群体
hp<-make_hp(hpsize=200,
            ng=500,h2=0.3,d2=0.2,phen_var=1,
            genome=genome,mutr=5*10**-4,sel_seq_qtl=0.1,sel_seq_mrk=0.05,laf=0.5)

## 第二步：从历史群体中抽样创建起始群体并进行扩充

# 首先是品种A，从历史群体中随机抽样公母各50头
Breed_A_Male_fndrs<-data.frame(number=50,select='rnd')
Breed_A_Female_fndrs<-data.frame(number=50,select='rnd')

# 从历史群体抽样之后通过连续随机交配来扩大群体，每次都是随机选择50头公畜与100头母畜来进行随机交配
Selection<-data.frame(matrix(NA, nrow=2, ncol=3))
names(Selection)<-c('Number','type','Value')
Selection$Number[1:2]<-c(50,100)
Selection$type[1:2]<-c('phen','phen')
Selection$Value[1:2]<-c('l','l')
Selection

# # 可以选择储存第N个世代的数据，这里选择了第10个世代
# Breed_A_data<-data.frame(matrix(NA, nrow=1, ncol=6))
# names(Breed_A_data)<-c('data','qtl','marker','seq','freq_qtl','freq_mrk')
# Breed_A_data[,1]<-10
# Breed_A_data[,2]<-10
# Breed_A_data[,3]<-10
# Breed_A_data[,4]<-10
# Breed_A_data[,5]<-10
# Breed_A_data[,6]<-10
# Breed_A_data

# litter_size设为5，从而达到扩充群体的目的
Breed_A<-sample_hp(hp_out=hp,Male_founders=Breed_A_Male_fndrs,
                   Female_founders=Breed_A_Female_fndrs,
                   ng=10,Selection=Selection,
                   litter_size=5,Display=FALSE)

# 接下来是品种B，也是从历史群体中随机抽样公母各50头
Breed_B_Male_fndrs<-data.frame(number=50,select='rnd')
Breed_B_Female_fndrs<-data.frame(number=50,select='rnd')

# 为了让两个品种的遗传结构的差异进一步加大，品种B按照表型进行从高到低的选择
Selection<-data.frame(matrix(NA, nrow=2, ncol=3))
names(Selection)<-c('Number','type','Value')
Selection$Number[1:2]<-c(50,100)
Selection$type[1:2]<-c('phen','phen')
Selection$Value[1:2]<-c('h','h') #这里表示优先挑选表型值较大的个体
Selection

# 品种B也是经过10个世代，每个世代每头母猪产生5头仔畜
Breed_B<-sample_hp(hp_out=hp,Male_founders=
                     Breed_B_Male_fndrs,Female_founders=Breed_B_Female_fndrs,
                   ng=10,Selection=Selection,
                   litter_size=5)

## 第三步：用于杂交的初始群体构建

# 从品种A中随机选择100头公畜、200头母畜，每头母畜产生5头仔畜，这样就构建了群体规模为1000头个体的A0
founder_pop1<-data.frame(matrix(NA, nrow=2, ncol=3))
names(founder_pop1)<-c('size','generation','select')
founder_pop1[1,]<-c(100,10,'rnd')
founder_pop1[2,]<-c(200,10,'rnd') 
founder_pop1

# 从品种B中随机选择100头公畜、200头，每头母畜产生5头仔畜，这样就构建了群体规模为1000头个体的B0
founder_pop2<-data.frame(matrix(NA, nrow=2, ncol=3))
names(founder_pop2)<-c('size','generation','select')
founder_pop2[1,]<-c(100,10,'rnd')
founder_pop2[2,]<-c(200,10,'rnd') 
founder_pop2

# 从品种A随机选择公畜100头，从品种B随机选择母畜200头，用来杂交产生AB0
founder_cross<-data.frame(matrix(NA, nrow=2, ncol=3))
names(founder_cross)<-c('pop','size','select')
founder_cross[1,]<-c('pop1',100,'rnd') # Select males from Breed A
founder_cross[2,]<-c('pop2',200,'rnd') # Select females from Breed B
founder_cross

## 第四步：杂交


# 针对A群体设置的选择模式，方法为GEBVC
Selection_pop1<-data.frame(matrix(NA, nrow=2, ncol=3))
names(Selection_pop1)<-c('Number','type','Value')
Selection_pop1$Number[1:2]<-c(100,200)
Selection_pop1$type[1:2]<-c('gebvc','gebvc')
Selection_pop1$Value[1:2]<-c('h','h')
Selection_pop1

# 同样地，针对B群体设置的选择方法也是GEBVC
Selection_pop2<-data.frame(matrix(NA, nrow=2, ncol=3))
names(Selection_pop2)<-c('Number','type','Value')
Selection_pop2$Number[1:2]<-c(100,200)
Selection_pop2$type[1:2]<-c('gebvc','gebvc')
Selection_pop2$Value[1:2]<-c('h','h')
Selection_pop2

# 设置如何从A、B两群体中选择杂交亲本的方式，方法依然为GEBVC
Cross_design<-data.frame(matrix(NA, nrow=2, ncol=4))
names(Cross_design)<-c('pop','size','select','value')
Cross_design[1,]<-c('pop1',100,'gebvc','h')
Cross_design[2,]<-c('pop2',200,'gebvc','h')
Cross_design

# 设置A群体的参考群体，大小为500头，随机选择（这里实际上相当于保留了全部个体500头）
train_A<-data.frame(matrix(NA, nrow=1, ncol=5))
names(train_A)<-c('size','sel','method','nIter','show')
train_A$size<-500
train_A$sel<-'rnd'
train_A$method<-'BayesB'
train_A$nIter<-1000 #抽样次数
train_A$show<-FALSE
train_A

# 针对B群体设置GEBV的估计方法
train_B<-data.frame(matrix(NA, nrow=1, ncol=5))
names(train_B)<-c('size','sel','method','nIter','show')
train_B$size<-500
train_B$sel<-'rnd'
train_B$method<-'BayesB'
train_B$nIter<-1000
train_B$show<-FALSE
train_B

# 保存AB1 - AB7的数据
output_cross<-data.frame(matrix(NA, nrow=7, ncol=5))
names(output_cross)<-c('data','qtl','marker','freq_qtl','freq_mrk')
output_cross[,1]<-1:7
output_cross[,2]<-1:7
output_cross[,3]<-1:7
output_cross[,4]<-1:7
output_cross[,5]<-1:7
output_cross

# 调取杂交函数进行杂交

set.seed(123)
cross_AB<-xbreed(pop1=Breed_A,pop2=Breed_B,founder_pop1=
                   founder_pop1,founder_pop2=founder_pop2,
                 founder_cross=founder_cross,
                 Selection_pop1=Selection_pop1,Selection_pop2=Selection_pop2,
                 Cross_design=Cross_design,train_type='purebred',
                 train_pop1=train_A,train_pop2=train_B,ng=7,litter_size=5,
                 saveAt='cross_pop_gebvc',output_cross=output_cross,Display=FALSE)

##如果我们设置GEBVP的选择模式来选择杂种的话，效果如何？

# 针对A群体设置的选择模式，方法为GEBVP
Selection_pop1<-data.frame(matrix(NA, nrow=2, ncol=3))
names(Selection_pop1)<-c('Number','type','Value')
Selection_pop1$Number[1:2]<-c(100,200)
Selection_pop1$type[1:2]<-c('gebv','gebv')
Selection_pop1$Value[1:2]<-c('h','h')
Selection_pop1

# 同样地，针对B群体设置的选择方法也是GEBVP
Selection_pop2<-data.frame(matrix(NA, nrow=2, ncol=3))
names(Selection_pop2)<-c('Number','type','Value')
Selection_pop2$Number[1:2]<-c(100,200)
Selection_pop2$type[1:2]<-c('gebv','gebv')
Selection_pop2$Value[1:2]<-c('h','h')
Selection_pop2

# 设置如何从A、B两群体中选择杂交亲本的方式，方法依然为GEBVC
Cross_design<-data.frame(matrix(NA, nrow=2, ncol=4))
names(Cross_design)<-c('pop','size','select','value')
Cross_design[1,]<-c('pop1',100,'gebv','h')
Cross_design[2,]<-c('pop2',200,'gebv','h')
Cross_design

set.seed(123)
cross_AB<-xbreed(pop1=Breed_A,pop2=Breed_B,founder_pop1=
                   founder_pop1,founder_pop2=founder_pop2,
                 founder_cross=founder_cross,
                 Selection_pop1=Selection_pop1,Selection_pop2=Selection_pop2,
                 Cross_design=Cross_design,train_type='purebred',
                 train_pop1=train_A,train_pop2=train_B,ng=7,litter_size=5,
                 saveAt='cross_pop_gebvp',output_cross=output_cross,Display=FALSE)

library(data.table)

dat1<-list()
for(i in 1:7){
  dat1[[i]]<-read.table(paste0("cross_pop_gebvc_cross__data_",i,".txt"),h=T)
}


dat2<-list()
for(i in 1:7){
  dat2[[i]]<-read.table(paste0("cross_pop_gebvp_cross__data_",i,".txt"),h=T)
}
sapply(dat2,function(x)mean(x$phen))

library(ggplot2)
library(ggpubr)
plot_dat<-data.frame(世代=rep(1:7,2),
                     杂种群体平均表型=c(sapply(dat1,function(x)mean(x$phen)),
                                sapply(dat2,function(x)mean(x$phen))),
                     方法=rep(c("GEBVC","GEBVP"),each=7))



png(file="结果.png",width=10, height=8,units='in',res=300)
ggplot(plot_dat, aes(x=世代, y=杂种群体平均表型,colour=方法)) + 
  geom_line() + geom_point() +
  theme_pubr(base_family="STSong", base_size = 15)

dev.off()

