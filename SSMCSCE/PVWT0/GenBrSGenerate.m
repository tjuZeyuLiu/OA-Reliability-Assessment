% 
BrS = [];
BrSNum = 0;
[BrS,ia,ic] = unique(GenBrS','rows','stable');
BrS = BrS';
BrSNum = size(BrS,2);
SPLoca = cell(BrSNum,1);
[B,numI] = sort(ic,'ascend');
i = 1;
istart = 1;
tic
while i<size(GenBrS,2)
    if B(i) == B(i+1)
        i = i + 1;
    else
        SPLoca{B(istart),1} = numI(istart:i);
        i = i + 1;
        istart = i;
    end       
end
SPLoca{B(istart),1} = numI(istart:i);
toc

savestr=strcat('114800SSMCSGenBrSdata20220905.mat');
save(savestr,'BrS','GenBrS','SPLoca','LoadS','BrSNum','W');