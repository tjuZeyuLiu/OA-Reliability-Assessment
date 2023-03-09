clear;
mpc0=load('case24');
mpc0=mpc0.mpc;

ldlv=load('ld17g33_lv8760');  %%从Lagrange的数据注意多一台发电机
ldlv=ldlv.ldlv;

 load('Ctgcase24.mat');

 load('pene0.05CEopt_ld17_cs24_lv8760_20220926_100000.mat');
CgOrigen = sparse(mpc0.origen(:,1),1:size(mpc0.origen,1),1,size(mpc0.gen,1),size(mpc0.origen,1));
GenBrU = CtgList{1}(:,3);
GenBrA = 1 - GenBrU;
NormalProb = prod(GenBrA);
optA = 1 - optU;


BrNum = size(mpc0.branch(:,1),1);
GenNum = size(mpc0.origen(:,1),1);
BusNum = size(mpc0.bus(:,1),1);
McsNum = 1000000;
ldlvnum = 8760;
datab  = zeros(58,ldlvnum);
tt = mpc0.bus~=0;
for i = 1: ldlvnum
    datab(tt,i) =mpc0.bus(tt) .* ldlv(i,1:17)';
     datab(BusNum*2+1:end,i) = CgOrigen *  (mpc0.origen(:,2).* ldlv(i,18:49)');
end
datab(BusNum+1:2*BusNum,:) = datab(1:BusNum,:);



% alphabet= 1:ldlvnum;
% prob = ldlv(:,50)./8760;
% LoadS = randsrc(McsNum,1,[alphabet;prob']);
% load('50000000SSMCS_ld17_cs24_lv100_20220408.mat');
LC=zeros(McsNum,1);
W=ones(McsNum,1);
beta =zeros(McsNum,1);
GenBrS = false(GenNum+BrNum,McsNum);
zeronum = 0;
datab0 = datab;
tic;
for i=1:McsNum

        mpc=mpc0;
        datab = datab0;
        GenBrRand = rand(GenNum+BrNum,1);
        GenBrS(:,i) = GenBrRand < optU;
        CtgListTmp = find(GenBrS(:,i) ==1);
        W(i,1) =  NormalProb ./ prod(GenBrA(CtgListTmp)) .* prod(GenBrU(CtgListTmp));
        optA = 1 - optU;
        W(i,1) =  W(i,1)./ (prod(optA) ./ prod(optA(CtgListTmp)) .* prod(optU(CtgListTmp)));
        LoadS(i,1) =  randi([1 8760]);
        ll = datab(:,LoadS(i,1));
       if (sum(GenBrS(:,i)) <= 1)
                  zeronum = zeronum +1;
        else
    
                    CtgCpntList = CpntList(GenBrS(:,i), :);
                     CtgGenList = (1 == CtgCpntList(:, 1));
                    CtgBrList = (2 == CtgCpntList(:, 1));
                     mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
                     mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
                     ll = datab(:,LoadS(i,1));
                     ll(2*BusNum+1:end,1)= CgOrigen * (mpc.origen(:,2) .* ldlv(LoadS(i,1),18:17+GenNum)');
%                      datab(2*BusNum+1:end,1:ldlvnum) = CgOrigen * (mpc.origen(:,2) .* ldlv(:,18:17+GenNum)');
                     mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
                    LC(i)=mcs_ld17_cal(mpc,ll);
      end
   if mod(i,1000000) == 0
       disp(i);
%        beta(i,1) = sqrt(var(LLC(1:i))/i)/mean(LLC(1:i))*100;
%        if beta(i,1) <= 1
%            disp('beta<1%');
%            break;
%        end
   end
end
time=toc
LLC = LC.*W; 
%   parfor i = 1:10000
% beta(i,1) = sqrt(var(LLC(1:i*100))/i/100)/mean(LLC(1:i*100))*100;
%    end

EENS = sum(LLC)/McsNum*8760*mpc0.baseMVA;
savestr=strcat('MCSCE_pene0.05_ld17_cs24_lv8760_20220926_1000000.mat');
save(savestr,'EENS','W','beta','LC','LLC','zeronum','time','GenBrS','LoadS');