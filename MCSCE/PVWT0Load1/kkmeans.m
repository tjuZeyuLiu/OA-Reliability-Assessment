function []=kkmeans(type,K)
str=strcat(type,'_lv8760.mat');

data=load(str);
data=data.ldlv;
data(:,end)=[];

X = data;
% X = [randn(50,2)+ones(50,2);randn(50,2)-ones(50,2);randn(50,2)+[ones(50,1),-ones(50,1)]];
opts = statset('Display','off'); 
%  […]=Kmeans(…,’Param1’,Val1,’Param2’,Val2,…)各输入输出参数介绍：     
%  X :N*P的数据矩阵       K: 表示将X划分为几类，为整数       Idx :N*1的向量，存储的是每个点的聚类标号       C: K*P的矩阵，存储的是K个聚类质心位置      sumD 1*K的和向量，存储的是类间所有点与该类质心点距离之和      D N*K的矩阵，存储的是每个点与所有质心的距离      […]=Kmeans(…,'Param1',Val1,'Param2',Val2,…)      这其中的参数Param1、Param2等，主要可以设置为如下：      1. ‘Distance’(距离测度)        ‘sqEuclidean’ 欧式距离（默认时，采用此距离方式）        ‘cityblock’ 绝度误差和，又称：L1        ‘cosine’ 针对向量        ‘correlation’  针对有时序关系的值        ‘Hamming’ 只针对二进制数据      2. ‘Start’（初始质心位置选择方法）        ‘sample’ 从X中随机选取K个质心点        ‘uniform’ 根据X的分布范围均匀的随机生成K个质心        ‘cluster’ 初始聚类阶段随机选择10%的X的子样本（此方法初始使用’sample’方法）         matrix 提供一K*P的矩阵，作为初始质心位置集合      3. ‘Replicates’（聚类重复次数）  整数

%调用Kmeans函数? ? ?X :N*P的数据矩阵 K: 表示将X划分为几类，为整数
%X N*P的数据矩阵
%Idx N*1的向量,存储的是每个点的聚类标号
%Ctrs K*P的矩阵,存储的是K个聚类质心位置
%SumD 1*K的和向量,存储的是类间所有点与该类质心点距离之和
%D N*K的矩阵，存储的是每个点与所有质心的距离; 
% K = 200;
% [Idx,Ctrs,SumD,D] = kmeans(X,50,'Replicates',100,'Options',opts); 
[Idx,Ctrs,SumD,D] = kmeans(X,K,'Replicates',200,'Options',opts); 

% [Idx,Ctrs] = kmeans(X,K); 
% load = zeros(
dim=size(X,2);
for i =  1 : K
    Ctrs(i,dim+1) = sum(Idx==i);
end
ldlv = Ctrs;

if dim == 1
ldlv=sortrows(ldlv,-1);
end
if dim == 17
    cc = [1.08;0.97;1.8;0.74;0.71;1.36;1.25;1.71;1.75;1.95;2.65;1.94;3.17;1;3.33;1.81;1.28];
    for i = 1 : K
          ldlv(i,dim+2)=cc'*ldlv(i,1:dim)';
    end
    ldlv=sortrows(ldlv,-(dim+2));
    ldlv(:,end) =[];
end

str=strcat(type,'_lv',num2str(K),'.mat');
save(str,'ldlv');
% legend('Cluster 1','Cluster 2','Cluster 3','Centroids','Location','NW');

% Ctrs

% SumD
% --------------------- 
% 作者：SethChai 
% 来源：CSDN 
% 原文：https://blog.csdn.net/a493823882/article/details/79282425 
% 版权声明：本文为博主原创文章，转载请附上博文链接！