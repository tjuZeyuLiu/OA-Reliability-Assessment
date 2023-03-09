clear;
clc;
warning off;
mpc0=load('case24');
mpc0=mpc0.mpc;
ldlv=load('ld17g33_lv8760');  %%从Lagrange的数据注意多一台发电机
ldlv=ldlv.ldlv;
ldlvnum = 8760;
 load('Ctgcase24.mat');
CgOrigen = sparse(mpc0.origen(:,1),1:size(mpc0.origen,1),1,size(mpc0.gen,1),size(mpc0.origen,1));
GenNum = size(mpc0.origen(:,1),1);
BusNum = size(mpc0.bus(:,1),1);
BrNum = size(mpc0.branch(:,1),1);

 datab  = zeros(58,ldlvnum);
 tt = mpc0.bus~=0;
for i = 1: ldlvnum
    datab(tt,i) =mpc0.bus(tt) .* ldlv(i,1:17)';
    datab(BusNum*2+1:end,i) = CgOrigen *  (mpc0.origen(:,2).* ldlv(i,18:49)');
end
datab(BusNum+1:2*BusNum,:) = datab(1:BusNum,:);


GenBrU = CtgList{1}(:,3);
GenBrA = 1 - GenBrU;



load('160200CEMCSGenBrSdata20220928.mat');
McsNum = 160200;
LC=zeros(BrSNum,1);
% beta =zeros(McsNum,1);
spnum =zeros(BrSNum,1);

zeronum = 0;
totsp = 50;
spnumdec =20;

%%%%% 提前建模
       mpc=mpc0;
       BusNum =size(mpc.bus,1);
       GenNum =size(mpc.gen,1);
       BrNum=size(mpc.branch,1);
       Ybus = sparse ([mpc.branch(:,1);mpc.branch(:,2)],[mpc.branch(:,2);mpc.branch(:,1)],[mpc.branch(:,3);mpc.branch(:,3)],BusNum,BusNum); %导纳矩阵
       Ybr = sparse ([1:BrNum,1:BrNum],[mpc.branch(:,1);mpc.branch(:,2)],[mpc.branch(:,3);-mpc.branch(:,3)],BrNum,BusNum);
       Yaa = - sum(Ybus);
       Ybus (sub2ind(size(Ybus),1:BusNum,1:BusNum)) = Yaa;
       Cg = sparse(mpc.gen(:,1),1:GenNum,1,BusNum,GenNum);
       EE = speye(2*BrNum+GenNum+BusNum);
       A = [ sparse(GenNum+BusNum,BusNum), speye(GenNum+BusNum);  
            Ybr ,                         sparse(BrNum,GenNum+BusNum);
            -Ybr,                         sparse(BrNum,GenNum+BusNum);];
        A = [A, EE];    
           A = [ Ybus , speye(BusNum), Cg, sparse(BusNum, 2 * BrNum + BusNum + GenNum);
             A;];
       b=[zeros(2*BusNum,1);mpc.gen(:,2);mpc.branch(:,4);mpc.branch(:,4)];
       c=[zeros(1,BusNum),ones(1,BusNum),zeros(1,BusNum + 2 * GenNum + 2 * BrNum)];
% %%%%% 提前建模
datab0 = datab;
databgen0 = datab(2*BusNum+1:end,1:ldlvnum);
databgen = databgen0;
datab(2*BusNum+1:end,:) = [];



tic;
for i=1:BrSNum

       mpc=mpc0;

        CtgCpntList = CpntList(GenBrS(:,numI(SPLoca(i,1))), :);
         CtgGenList = (1 == CtgCpntList(:, 1));
        CtgBrList = (2 == CtgCpntList(:, 1));
         FailGenList = CtgCpntList(CtgGenList, 2);

         if (SPLoca(i,2) - SPLoca(i,1))>100
                databgen = databgen0 - CgOrigen(:,FailGenList) * (mpc.origen(FailGenList,2) .* ldlv(:,17+FailGenList)');
         else
                LoadStt = LoadS(numI(SPLoca(i,1):SPLoca(i,2)));
                databgen(:,LoadStt) = databgen0(:,LoadStt) - CgOrigen(:,FailGenList) * (mpc.origen(FailGenList,2) .* ldlv(LoadStt,17+FailGenList)');
         end
          FailBrList = CtgCpntList(CtgBrList, 2);

         if  i == 1
                [LC(i,1),spnum(i,1),sp]=ssmcs_ld17_cal_init(A,b,c,FailBrList,mpc,numI(SPLoca(i,1):SPLoca(i,2)),LoadS,W,datab,databgen,totsp,CtgCpntList(CtgBrList, 2),spnumdec);
         else
                [LC(i,1),spnum(i,1),sp]=ssmcs_ld17_cal(A,b,c,FailBrList,mpc,numI(SPLoca(i,1):SPLoca(i,2)),LoadS,W,datab,databgen,sp,spnum(i-1),totsp,CtgCpntList(CtgBrList, 2),spnumdec);
         end


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
savestr=strcat('SSCEMCS_ld17_cs24_lv8760_160200_20220928.mat');
save(savestr,'LC','zeronum','time','GenBrS','LoadS','EENS','spnum');