% function []=sp_ld17_main(casenum,ldlvnum)
clear;
clc;
casenum = 24;
casenum=num2str(casenum);
casestr=strcat('case',casenum);
mpc0=load(casestr);
mpc0=mpc0.mpc;
ldlvnum = 100;
% ldlvnum=num2str(ldlvnum);
% ldlvstr=strcat('ld17_lv',ldlvnum);
% ldlv=load(ldlvstr);
% ldlv=ldlv.ldlv;
% ldlvnum = 100;

CtgLevelMax = 5;
load('Ctgcase24.mat');
load('pene0.05caserts79datablv100.mat');
% [ CtgList, CpntList ] = CreatSECtgList(mpc0, CtgLevelMax);



CgOrigen = sparse(mpc0.origen(:,1),1:size(mpc0.origen,1),1,size(mpc0.gen,1),size(mpc0.origen,1));
GenNum = size(mpc0.origen(:,1),1);
BusNum = size(mpc0.bus(:,1),1);
datab0 = datab;
tic;
[lc0,spnum0]=nor_ld17g33_cal(mpc0,ldlv,datab);

CtgListTmp = CtgList{1};
CtgNum = size(CtgListTmp,1);
lc1=zeros(CtgNum ,1);
spnum1=zeros(CtgNum ,1);
for i=1:CtgNum  
   mpc=mpc0;
   datab = datab0;
    CtgCpntNo = CtgListTmp(i,1);       % 故障支路详细列表
     CtgCpntList = CpntList(CtgCpntNo, :);
      CtgGenList = (1 == CtgCpntList(:, 1));
      CtgBrList = (2 == CtgCpntList(:, 1));
        mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
        mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
        datab(2*BusNum+1:end,1:ldlvnum) = CgOrigen * (mpc.origen(:,2) .* ldlv(:,18:17+GenNum)');
        mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
       [lc1(i),spnum1(i)]=nor_ld17g33_cal(mpc,ldlv,datab);
end

CtgListTmp = CtgList{2};
CtgNum = size(CtgListTmp,1);
lc2=zeros(CtgNum ,1);
spnum2=zeros(CtgNum ,1);
for i=1:CtgNum  %% 
   mpc=mpc0;
    datab = datab0;
    CtgCpntNo = CtgListTmp(i,1:2);       % 故障支路详细列表
     CtgCpntList = CpntList(CtgCpntNo, :);
      CtgGenList = (1 == CtgCpntList(:, 1));
      CtgBrList = (2 == CtgCpntList(:, 1));
        mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
        mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
        datab(2*BusNum+1:end,1:ldlvnum) = CgOrigen * (mpc.origen(:,2) .* ldlv(:,18:17+GenNum)');
        mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
       [lc2(i),spnum2(i)]=nor_ld17g33_cal(mpc,ldlv,datab);
end

CtgListTmp = CtgList{3};
CtgNum = size(CtgListTmp,1);
lc3=zeros(CtgNum ,1);
spnum3=zeros(CtgNum ,1);
for i=1:CtgNum  %% 
   mpc=mpc0;
   datab = datab0;
    CtgCpntNo = CtgListTmp(i,1:3);       % 故障支路详细列表
     CtgCpntList = CpntList(CtgCpntNo, :);
      CtgGenList = (1 == CtgCpntList(:, 1));
      CtgBrList = (2 == CtgCpntList(:, 1));
        mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
        mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
        mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
          datab(2*BusNum+1:end,1:ldlvnum) = CgOrigen * (mpc.origen(:,2) .* ldlv(:,18:17+GenNum)');
       [lc3(i),spnum3(i)]=nor_ld17g33_cal(mpc,ldlv,datab);
end

CtgListTmp = CtgList{4};
CtgNum = size(CtgListTmp,1);
lc4=zeros(CtgNum ,1);
spnum4=zeros(CtgNum ,1);
for i=1:CtgNum  %% 
   mpc=mpc0;
   datab = datab0;
    CtgCpntNo = CtgListTmp(i,1:4);       % 故障支路详细列表
     CtgCpntList = CpntList(CtgCpntNo, :);
      CtgGenList = (1 == CtgCpntList(:, 1));
      CtgBrList = (2 == CtgCpntList(:, 1));
        mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
        mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
         datab(2*BusNum+1:end,1:ldlvnum) = CgOrigen * (mpc.origen(:,2) .* ldlv(:,18:17+GenNum)');
        mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
        
       [lc4(i),spnum4(i)]=nor_ld17g33_cal(mpc,ldlv,datab);
end

CtgListTmp = CtgList{5};
CtgNum = size(CtgListTmp,1);
lc5 = zeros(CtgNum ,1);
spnum5 = zeros(CtgNum ,1);
for i=1:CtgNum  %%
   mpc=mpc0;
      datab = datab0;
    CtgCpntNo = CtgListTmp(i,1:5);       % 故障支路详细列表
     CtgCpntList = CpntList(CtgCpntNo, :);
      CtgGenList = (1 == CtgCpntList(:, 1));
      CtgBrList = (2 == CtgCpntList(:, 1));
        mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
        mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
         datab(2*BusNum+1:end,1:ldlvnum) = CgOrigen * (mpc.origen(:,2) .* ldlv(:,18:17+GenNum)');
        mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
       [lc5(i),spnum5(i)]=nor_ld17g33_cal(mpc,ldlv,datab);
end

time=toc
sumspnum=spnum0 + sum(spnum1) + sum(spnum2) + sum(spnum3) + sum(spnum4) + sum(spnum5);
lc{1}=lc1;
lc{2}=lc2;
lc{3}=lc3;
lc{4}=lc4;
lc{5}=lc5;
SECalcullateReliabilityIndices;
% spnum=mean([spnum0;spnum1;spnum2;spnum3;spnum4;spnum5]);

savestr=strcat('nor_ld17g33pene0.05_cs',casenum,'_lv',ldlvnum,'20220927','.mat');
save(savestr,'lc0','lc','LC','LC1','LC2','LC3','LC4','LC5','IISELC','IISELC2','IISELC3','IISELC4','IISELC5','sumspnum','spnum0','spnum1','spnum2','spnum3','spnum4','spnum5','time');
