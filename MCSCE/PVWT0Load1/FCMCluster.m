
function[center,U,obj_fun]=FCMCluster(data,n,options)

%采用模糊C均值将数据集data分为n类

%用法

% 1 [center,U,obj_fcn]=FCMCluster(data,n,options);

% 2 [center,U,obj_fcn]=FCMCluster(data,n);

 

%输入 

% data    n*m矩阵，n个样本数，每个样本的维度为m

% n       类别数

% options 4*1 矩阵

%   options(1):隶属度矩阵U的加权指数

%   options(2):最大迭代次数

%   options(3):隶属度最小变化量，迭代终止条件

%   options(4):每次迭代是否输出信息标志

 

%输出

% center    聚类中心

% U         隶属度矩阵

% obj_fun   目标函数值

 

 

if nargin~=2 && nargin~=3

    error('Too many or too few input arguments');

end 

 

data_n=size(data,1);

in_n=size(data,2);

 

%默认参数

default_options=[3;1000;1e-4;1];

 

%参数配置

  %如果只输入前两个参数，选用默认的参数;如果参数个数小于4，其他选用默认参数

  if nargin==2

      options=default_options;

  else

       if length(options)<4

           tmp=default_options;

           tmp(1:length(options))=options;

           options=tmp;

       end 

       nan_index=find(isnan(options)==1);

       options(nan_index)=default_options(nan_index);

 

       if options(1)<=1

           error('The exponent should be greater than 1!');

       end 

  end 

 

  %将options 中的分量分别赋值给四个变量

  expo=options(1);

  max_iter=options(2);

  min_impro=options(3);

  display=options(4);

 

  obj_fun=zeros(max_iter,1);

 

  %初始化模糊分配矩阵

  U=initfcm(n,data_n);

 

 

  %主程序

   for i=1:max_iter

       [U,center,obj_fun(i)]=stepfcm(data,U,n,expo);

       if display

           fprintf('FCM:Iteration count=%d,obj_fun=%f\n',i,obj_fun(i));

       end

       %终止条件判别

       if i>1

           if abs(obj_fun(i)-obj_fun(i-1))<min_impro

               break;

           end 

       end 

   end

   iter_n=i;

   obj_fun(iter_n+1:max_iter)=[];

%   

 

end

%%子函数 模糊矩阵初始化

    function U= initfcm(n,data_n)

        U=rand(n,data_n);

        col_sum=sum(U);

        U=U./col_sum(ones(n,1),:);

    end

 

 

%%子函数 逐步聚类

        function [U_new,center,obj_fun]=stepfcm(data,U,n,expo)

            mf=U.^expo;

            center=mf*data./((ones(size(data,2),1)*sum(mf'))');

            dist=distfcm(center,data);

            obj_fun=sum(sum((dist.^2).*mf));

            tmp=dist.^(-2/(expo-1));

            U_new=tmp./(ones(n,1)*sum(tmp));

        end

 

 

%%子函数 计算距离

            function out=distfcm(center,data)

                out=zeros(size(center,1),size(data,1));

                for k=1:size(center,1)

                    out(k,:)=sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));





                end 

            end
