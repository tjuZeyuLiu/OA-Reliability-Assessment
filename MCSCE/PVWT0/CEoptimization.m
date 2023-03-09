clear;
mpc0=load('case24');
mpc0=mpc0.mpc;

ldlv=load('ld17_lv8760');
ldlv=ldlv.ldlv;
ldlvnum = 8760;
 load('Ctgcase24.mat');
%  load('case24datablv8760.mat');
 
CgOrigen = sparse(mpc0.origen(:,1),1:size(mpc0.origen,1),1,size(mpc0.gen,1),size(mpc0.origen,1));
GenBrU = CtgList{1}(:,3);
GenBrA = 1 - GenBrU;
NormalProb = prod(GenBrA);
BrNum = size(mpc0.branch(:,1),1);
GenNum = size(mpc0.origen(:,1),1);
BusNum = size(mpc0.bus(:,1),1);
datab  = zeros(48,ldlvnum);
 tt = mpc0.bus~=0;
for i = 1: ldlvnum
    datab(tt,i) =mpc0.bus(tt) .* ldlv(i,1:17)';
end
datab(BusNum+1:2*BusNum,:) = datab(1:BusNum,:);


OptMaxNum = 10000;
McsNum = 100000;
McsNp = 10000;

% alphabet= 1:ldlvnum;
% prob = ldlv(:,18)./8760;
% LoadS = randsrc(McsNum,1,[alphabet;prob']);
% load('50000000SSMCS_ld17_cs24_lv100_20220408.mat');
LC=zeros(McsNum,1);
beta =zeros(McsNum,1);
LoadS =zeros(McsNum,1);
GenBrS = false(GenNum+BrNum,McsNum);
zeronum = 0;
tic;
parfor i=1:McsNum

        mpc=mpc0;
        GenBrRand = rand(GenNum+BrNum,1);
        GenBrS(:,i) = GenBrRand < GenBrU;
        LoadS(i,1) =  randi([1 8760]);
        ll = datab(:,LoadS(i,1));
        if (sum(GenBrS(:,i)) <= 1)
                  zeronum = zeronum +1;
                   CtgCpntList = CpntList(GenBrS(:,i), :);
                     CtgGenList = (1 == CtgCpntList(:, 1));
                    CtgBrList = (2 == CtgCpntList(:, 1));
                     mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
                     mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
                     mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
                    LC(i) = sum(mpc.gen(:,2));
        else
                
                    CtgCpntList = CpntList(GenBrS(:,i), :);
                     CtgGenList = (1 == CtgCpntList(:, 1));
                    CtgBrList = (2 == CtgCpntList(:, 1));
                     mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
                     mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
                     mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
                    LC(i)=mcs_ld17_cal(mpc,ll);
                   LC(i) = sum(mpc.gen(:,2));
      end
   if mod(i,1000000) == 0
       disp(i);
   end
end

rk =zeros(OptMaxNum,1);
optU = zeros(GenNum+BrNum,OptMaxNum);
W=ones(McsNum,1);

k = 1;
optU(:,k) = GenBrU;

[S,IX] = sort(LC(1:McsNum));
LOLP=zeros(McsNum,1);
flag = 0;
if S(McsNp)>28.5
% if S(McsNp)<0.01
    rk(k,1) = flag;
    LOLP(IX(1:McsNp))= 1;
else
    flag = 1;
    rk(k,1) = flag;
    LOLP = LC(1:McsNum) <28.5;
%       LOLP = LC(1:McsNum) >0.00001;
end
for i=1:McsNum
    CtgListTmp = find(GenBrS(:,i) ==1);
    W(i,1) =  NormalProb ./ prod(GenBrA(CtgListTmp)) .* prod(GenBrU(CtgListTmp));
    optA = 1 - optU(:,k);
    W(i,1) =  W(i,1)./ (prod(optA) ./ prod(optA(CtgListTmp)) .* prod(optU(CtgListTmp,k)));
end
k = k + 1;
optU(:,k) = 0.9999* (1-sum(LOLP'.*W'.*(1-GenBrS(:,1:McsNum)),2)./sum(LOLP.*W));
while flag~= 1
    parfor i=1:McsNum

        mpc=mpc0;
        GenBrRand = rand(GenNum+BrNum,1);
        GenBrS(:,i) = GenBrRand < optU(:,k);
        LoadS(i,1) =  randi([1 8760]);
        ll = datab(:,LoadS(i,1));
        if (sum(GenBrS(:,i)) <= 1)
                  zeronum = zeronum +1;
                     CtgCpntList = CpntList(GenBrS(:,i), :);
                     CtgGenList = (1 == CtgCpntList(:, 1));
                    CtgBrList = (2 == CtgCpntList(:, 1));
                     mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
                     mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
                     mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
                    LC(i) = sum(mpc.gen(:,2));
        else
                
                    CtgCpntList = CpntList(GenBrS(:,i), :);
                     CtgGenList = (1 == CtgCpntList(:, 1));
                    CtgBrList = (2 == CtgCpntList(:, 1));
                     mpc.origen(CtgCpntList(CtgGenList, 2),2) = 0;
                     mpc.gen(:,2) = CgOrigen *  mpc.origen(:,2);
                     mpc.branch(CtgCpntList(CtgBrList, 2), [3,4]) = 0;
                    LC(i)=mcs_ld17_cal(mpc,ll);
                     LC(i) = sum(mpc.gen(:,2));
        end
    end
    
    [S,IX] = sort(LC(1:McsNum));
    LOLP=zeros(McsNum,1);
    flag = 0;
if S(McsNp)>28.5
% if S(McsNp)<0.01
    rk(k,1) = flag;
    LOLP(IX(1:McsNp))= 1;
else
    flag = 1;
    rk(k,1) = flag;
    LOLP = LC(1:McsNum) <28.5;
%       LOLP = LC(1:McsNum) >0.00001;
end
    for i=1:McsNum
        CtgListTmp = find(GenBrS(:,i) ==1);
        W(i,1) =  NormalProb ./ prod(GenBrA(CtgListTmp)) .* prod(GenBrU(CtgListTmp));
        optA = 1 - optU(:,k);
        W(i,1) =  W(i,1)./ (prod(optA) ./ prod(optA(CtgListTmp)) .* prod(optU(CtgListTmp,k)));
    end
    k = k + 1;
    optU(:,k) = 0.9999* (1-sum(LOLP'.*W'.*(1-GenBrS(:,1:McsNum)),2)./sum(LOLP.*W))+0.0001*optU(:,k-1);
    
    disp(k);
    
end

time=toc

optU = optU(:,k);
% EENS = sum(LC)/McsNum*8760*mpc0.baseMVA;
savestr=strcat('CEopt_ld17_cs24_lv8760_20220926_100000.mat');
save(savestr,'optU','time');