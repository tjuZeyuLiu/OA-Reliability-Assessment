function [plc,sp]=sp_newmskopt(a,b,c,BusNum,GenNum,CtgBrList)
prob.a=a;
prob.c=c';
prob.blc=b;
prob.buc=b;
prob.blx=zeros(size(c,2),1);
prob.bux=[];
cmd='minimize echo(0)';
mosek_opt.MSK_IPAR_NUM_THREADS = 1;
% [x,fval,exitflag,output,lambda] = linprog(c',[],[],a,b,zeros(size(c,2),1),[],[],[]);
% mosek_opt.MSK_IPAR_INTPNT_MULTI_THREAD = 'MSK_ON';
[~,res]=mosekopt(cmd,prob,mosek_opt);
% [~,res]=mosekopt(cmd,prob);
plc=res.sol.bas.pobjval;
sp.xbflag = res.sol.bas.skx(:,1)=='B';
cb=find(res.sol.bas.skc(:,1)=='B');


if ~isempty(cb)
   
   sp.xbflag(cb(cb<=BusNum)+BusNum) = 1;  %% 形成孤岛，无负荷节点，RTS79中节点24代表，对应负荷削减为0进入基，不然只有松弛入基容易奇异
   sp.xbflag(cb(cb>2*BusNum+GenNum)+BusNum+GenNum) = 1;  
end
sp.A = a;
sp.xb = find(sp.xbflag == 1);
sp.xn = find(sp.xbflag == 0);

sp.judge=find((sp.xb>BusNum&sp.xb<=2*BusNum)|sp.xb>2*BusNum+GenNum);
sp.B = a(:,sp.xb);
sp.invB=inv(sp.B);
sp.judgeinvB=sp.invB(sp.judge,:);
% sp.invB=invB(BusNum+1:end,:); %% error 所求x为最优基解，busnum部分很可能本来就没有
sp.w=c(sp.xb)*sp.invB;
sp.CtgBrList = CtgBrList;
sp.CtgBrNum = size(CtgBrList,1);     
sp.num = 1;