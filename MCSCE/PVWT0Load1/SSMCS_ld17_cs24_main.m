clear;
clc;
warning off;
mpc0=load('case24');
mpc0=mpc0.mpc;

ldlv=load('ld17_lv100');
ldlv=ldlv.ldlv;

 load('Ctgcase24.mat');
 load('case24datablv100.mat');
CgOrigen = sparse(mpc0.origen(:,1),1:size(mpc0.origen,1),1,size(mpc0.gen,1),size(mpc0.origen,1));
GenBrU = CtgList{1}(:,3);
GenBrA = 1 - GenBrU;
BrNum = size(mpc0.branch(:,1),1);
GenNum = size(mpc0.origen(:,1),1);
McsNum = 25000000;
ldlvnum = 100;
% alphabet= 1:ldlvnum;
% prob = ldlv(:,18)./8760;
% LoadS = randsrc(McsNum,1,[alphabet;prob']);
load('25000000SSMCSGenBrSdata20220514.mat');
% % load('spnum-5SSMCS_ld17_cs24_lv100_20220408.mat');
LC=zeros(BrSNum,1);
% beta =zeros(McsNum,1);
spnum =zeros(BrSNum,1);
% GenBrS = false(GenNum+BrNum,McsNum);
zeronum = 0;

totsp = 8;
tic;
for i=1:BrSNum

        mpc=mpc0;

        CtgCpntList = CpntList(GenBrS(:,SPLoca{i}(1)), :);
         CtgGenList = (1 == CtgCpntList(:, 1));
        CtgBrList = (2 == CtgCpntList(:, 1));
         mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
         mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
         mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
         if  i == 1
                [LC(i,1),spnum(i,1),sp]=ssmcs_ld17_cal_init(mpc,SPLoca{i},LoadS,datab,totsp,CtgCpntList(CtgBrList, 2));
         else
                [LC(i,1),spnum(i,1),sp]=ssmcs_ld17_cal(mpc,SPLoca{i},LoadS,datab,sp,spnum(i-1),totsp,CtgCpntList(CtgBrList, 2));
         end
%       if mod(i,10) == 0
% for j = 1 : totsp
%     sp(j).num =0;
% end
%          end

   if mod(i,1000000) == 0
       disp(i);
%        beta(i,1) = sqrt(var(LC(1:i))/i)/mean(LC(1:i))*100;
%        if beta(i,1) <= 1
%            disp('beta<1%');
%            break;
%        end
   end
end

time=toc
EENS = sum(LC)/McsNum*8760*mpc0.baseMVA;
savestr=strcat('spnum-4SSMCS_ld17_cs24_lv100_25000000_20220524.mat');
save(savestr,'LC','zeronum','time','GenBrS','LoadS','EENS','spnum');