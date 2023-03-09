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
load('pene0.15caserts79datablv100.mat');
% [ CtgList, CpntList ] = CreatSECtgList(mpc0, CtgLevelMax);

CtgNum = size(CtgList{1},1);
lc{1}=zeros(CtgNum ,1);

CtgNum=size(CtgList{2},1);
lc{2}=zeros(CtgNum ,1);

CtgNum=size(CtgList{3},1);
lc{3}=zeros(CtgNum ,1);

CtgNum=size(CtgList{4},1);
lc{4}=zeros(CtgNum ,1);

CtgNum=size(CtgList{5},1);
lc{5}=zeros(CtgNum ,1);

GenNum = size(mpc0.origen(:,1),1);
BusNum = size(mpc0.bus(:,1),1);
CgOrigen = sparse(mpc0.origen(:,1),1:size(mpc0.origen,1),1,size(mpc0.gen,1),size(mpc0.origen,1));
totsp = 20;
spnum1=zeros(size(CtgBr{1},1) ,1);
spnum2=zeros(size(CtgBr{2},1) ,1);
spnum3=zeros(size(CtgBr{3},1) ,1);
datab0 =datab;
tic;
[lc0,spnum0,sp]=sp_ld17g33_cal_init(mpc0,ldlv,datab,totsp,[]);

% % 0阶线路故障

for i=1:size(CtgBr{1},1)
        mpc=mpc0;
        datab = datab0;
        CtgCpntNo = CtgList{CtgBr{1}(i,1)}(CtgBr{1}(i,2),1:CtgBr{1}(i,1));   % 故障支路详细列表
        CtgCpntList = CpntList(CtgCpntNo, :);
        CtgGenList = (1 == CtgCpntList(:, 1));
        CtgBrList = (2 == CtgCpntList(:, 1));
        mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
        mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
          datab(2*BusNum+1:end,1:ldlvnum) = CgOrigen * (mpc.origen(:,2) .* ldlv(:,18:17+GenNum)');
        mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
        if  i == 1
            [lc{CtgBr{1}(i,1)}(CtgBr{1}(i,2)),spnum1(i),sp]=sp_ld17g33_cal(mpc,ldlv,datab,sp,spnum0,totsp,CtgCpntList(CtgBrList, 2));
        else
           [lc{CtgBr{1}(i,1)}(CtgBr{1}(i,2)),spnum1(i),sp]=sp_ld17g33_cal(mpc,ldlv,datab,sp,spnum1(i-1),totsp,CtgCpntList(CtgBrList, 2));
        end
        
     if (mod(i,30000)==0)
       disp(i);
   end
end
sumspnum1 = spnum1(size(CtgBr{1},1));
% toc
% 
% % % % 1阶线路故障
sumspnum2 = 0;
for i= 1:1:size(CtgBr{2},1)
    mpc=mpc0;
     datab = datab0;
    CtgCpntNo = CtgList{CtgBr{2}(i,1)}(CtgBr{2}(i,2),1:CtgBr{2}(i,1));   % 故障支路详细列表
    CtgCpntList = CpntList(CtgCpntNo, :);
    CtgGenList = (1 == CtgCpntList(:, 1));
    CtgBrList = (2 == CtgCpntList(:, 1));
    mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
    mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
    datab(2*BusNum+1:end,1:ldlvnum) = CgOrigen * (mpc.origen(:,2) .* ldlv(:,18:17+GenNum)');
    mpc.branch(CtgCpntList(CtgBrList, 2),[3,4]) = 0;
    if i == 1
           [lc{CtgBr{2}(i,1)}(CtgBr{2}(i,2)),spnum2(i),sp]=sp_ld17g33_cal(mpc,ldlv,datab,sp,spnum1(size(CtgBr{1},1)),totsp,CtgCpntList(CtgBrList, 2));
    else
           [lc{CtgBr{2}(i,1)}(CtgBr{2}(i,2)),spnum2(i),sp]=sp_ld17g33_cal(mpc,ldlv,datab,sp,spnum2(i-1),totsp,CtgCpntList(CtgBrList, 2));
    end
   if (mod(i,30000)==0)
       disp(i);
   end
end
sumspnum2 = spnum2(size(CtgBr{2},1))-sumspnum1;

% % 2阶线路故障
for i= 1:size(CtgBr{3},1)
        mpc=mpc0;
         datab = datab0;
        CtgCpntNo = CtgList{CtgBr{3}(i,1)}(CtgBr{3}(i,2),1:CtgBr{3}(i,1));   % 故障支路详细列表
        CtgCpntList = CpntList(CtgCpntNo, :);
        CtgGenList = (1 == CtgCpntList(:, 1));
        CtgBrList = (2 == CtgCpntList(:, 1));
        mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
        mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
        datab(2*BusNum+1:end,1:ldlvnum) = CgOrigen * (mpc.origen(:,2) .* ldlv(:,18:17+GenNum)');
       mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
     if i == 1
           [lc{CtgBr{3}(i,1)}(CtgBr{3}(i,2)),spnum3(i),sp]=sp_ld17g33_cal(mpc,ldlv,datab,sp,spnum2(size(CtgBr{2},1)),totsp,CtgCpntList(CtgBrList, 2));
    else
           [lc{CtgBr{3}(i,1)}(CtgBr{3}(i,2)),spnum3(i),sp]=sp_ld17g33_cal(mpc,ldlv,datab,sp,spnum3(i-1),totsp,CtgCpntList(CtgBrList, 2));
    end
   if (mod(i,30000)==0)
       disp(i);
   end
end
sumspnum3 = spnum3(size(CtgBr{3},1))-spnum2(size(CtgBr{2},1));

time=toc
sumspnum=spnum0 + sumspnum1 + sumspnum2 + sumspnum3;
SECalcullateReliabilityIndices;
% spnum=mean([spnum0;spnum1;spnum2;spnum3;spnum4;spnum5]);

savestr=strcat('pene0.15sp_ld17_cs',casenum,'_lv',ldlvnum,'20220928','.mat');
save(savestr,'sp','lc0','lc','LC','LC1','LC2','LC3','LC4','LC5','IISELC','IISELC2','IISELC3','IISELC4','IISELC5','sumspnum','spnum0','spnum1','spnum2','spnum3','sumspnum1','sumspnum2','sumspnum3','time');
