function plc=nor_mskopt(a,b,c)
prob.a=a;
prob.c=c';
prob.blc=b;
prob.buc=b;
prob.blx=zeros(size(c,2),1);
prob.bux=[];
cmd='minimize echo(0)';
mosek_opt.MSK_IPAR_NUM_THREADS = 1;
[~,res]=mosekopt(cmd,prob,mosek_opt);
plc=res.sol.bas.pobjval;
% options = optimoptions('linprog','Algorithm','interior-point','Display','off');
% lb = zeros(size(c,2),1);
% [x,fval] = linprog(c,[],[],a,b,lb,[],options);
% plc=fval;