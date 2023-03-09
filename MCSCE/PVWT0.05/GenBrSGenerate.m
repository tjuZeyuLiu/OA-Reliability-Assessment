% 
GenBrS = GenBrS(:,1:185300);
W = W(1:185300);
LoadS = LoadS(1:185300);
BrS = [];
BrSNum = 0;
[BrS,ia,ic] = unique(GenBrS','rows','stable');
BrS = GenBrS(:,ia);
% BrS = BrS';
BrSNum = size(BrS,2);
SPLoca = cell(BrSNum,1);
[B,numI] = sort(ic,'ascend');
i = 1;
istart = 1;
SPLoca = zeros(BrSNum,2);
[BB,ia,ic] = unique(B,'rows','stable');
SPLoca(:,1) = ia;
SPLoca(1:BrSNum-1,2) = ia(2:BrSNum)-1;
SPLoca(BrSNum,2) = size(GenBrS,2);

savestr=strcat('185300CEMCSGenBrSdata20220926.mat');
save(savestr,'numI','BrS','GenBrS','SPLoca','LoadS','BrSNum','W');