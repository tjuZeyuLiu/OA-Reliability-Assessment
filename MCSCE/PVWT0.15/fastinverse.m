function [invB]=fastinverse1(reinvB,reB,B,BNum)

deltaB = B - reB;

row = find(sum(abs(deltaB),2)~=0);  %行求和
nrow = size(row,1);
BB = sparse(row,1:nrow,ones(nrow,1),BNum,nrow,nrow);
CC = deltaB(row,:);
cc = (speye(nrow,nrow)+CC*reinvB*BB)\CC;
invB = reinvB - reinvB*BB*cc*reinvB;