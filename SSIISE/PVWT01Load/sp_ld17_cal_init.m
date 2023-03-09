function [lc,spnum,sp]=sp_ld17_cal_init(mpc,ldlv,datab,totsp,CtgBrList)
   BusNum =size(mpc.bus,1);
   GenNum =size(mpc.gen,1);
   BrNum=size(mpc.branch,1);
   ldlvnum = size(ldlv,1);
%    B=zeros(BrNum,1);
%    A=sparse(2*BusNum + GenNum + 2*BrNum,3*BusNum + 2*GenNum + 2*BrNum);
%    A= [bus     bus      bus    gen   gen  branch branch
%     bus
%     gen
%     branch 
%     branch] 
   Ybus = sparse ([mpc.branch(:,1);mpc.branch(:,2)],[mpc.branch(:,2);mpc.branch(:,1)],[mpc.branch(:,3);mpc.branch(:,3)],BusNum,BusNum); %导纳矩阵
   Ybr = sparse ([1:BrNum,1:BrNum],[mpc.branch(:,1);mpc.branch(:,2)],[mpc.branch(:,3);-mpc.branch(:,3)],BrNum,BusNum);
   Yaa = - sum(Ybus);
% %    Yaa(Yaa==0) = 1;
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
%       if ~isempty(CtgBrList)
%        b(BusNum*2+GenNum+ BrNum+CtgBrList,:) = [];
%        b(BusNum*2+GenNum + CtgBrList,:) = [];
%        A(BusNum*2+GenNum+ BrNum + CtgBrList,:) = [];
%        A(BusNum*2+GenNum+CtgBrList,:) = [];
%        A(:,BusNum*3+GenNum*2+ BrNum + CtgBrList) = [];
%        A(:,BusNum*3+GenNum*2+CtgBrList) = [];
%        BrNum = BrNum - size(CtgBrList,1);
%    end
   c=[zeros(1,BusNum),ones(1,BusNum),zeros(1,BusNum + 2 * GenNum + 2 * BrNum)];
   lc=0;

   sp(1:totsp) = struct('invB',[],'w',[],'xb',[],'xn',[],'CtgBrList',[],'CtgBrNum',[],'judge',[],'xbflag',[],'num',[]);
   spnum=0;
   CtgBrNum = size(CtgBrList,1);
   for l=1:ldlvnum
      b(1:2*BusNum)= datab(:,l);

% datab = zeros(48,100);
% for i = 1 : ldlvnum
%       datab(mpc.area,i)=mpc.bus(mpc.area).*ldlv(i,1:17)';
%       datab(BusNum+1:2*BusNum,i)=datab(1:BusNum,i);
% end
      flag=0;
      for i=spnum:-1:max(spnum-1,1)
          if ((sp(mod(i-1,totsp)+1).CtgBrNum == 0)&& (CtgBrNum == 0))|(sp(mod(i-1,totsp)+1).CtgBrList == CtgBrList)
          if sp(mod(i-1,totsp)+1).invB*b > -1e-8
             plc=sp(mod(i-1,totsp)+1).w*b;
             sp(mod(i-1,totsp)+1).num = sp(mod(i-1,totsp)+1).num + 1;
             flag=1;
             break;
          end
          end
      end
      if flag==0
         spnum=spnum+1;
         [plc,sp(mod(spnum-1,totsp)+1)]=sp_newmskopt(A,b,c,BusNum,GenNum,CtgBrList);
    
      end
      lc = lc + plc*ldlv(l,2);
end