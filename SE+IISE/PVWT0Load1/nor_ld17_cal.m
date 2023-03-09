function [lc,spnum,sp]=nor_ld17_cal(mpc,ldlv,datab,sp,spnum,totsp,CtgBrList)
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
%    islands = find(Yaa == 0);
%    Yaa(islands) = 1;
   Ybus (sub2ind(size(Ybus),1:BusNum,1:BusNum)) = Yaa;
%    if ~isempty(islands)
%        for i = 1 : size(islands,1)
%                CtgBrislands = find( (mpc.branch(CtgBrList,1) == islands(i,1)) |  (mpc.branch(CtgBrList,2) == islands(i,1)) );
%                Ybr(CtgBrList(CtgBrislands(1)),islands(i)) = 1;
%        end
%    end
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
%        BrNum = BrNum -  size(CtgBrList,1);
%    end
   c=[zeros(1,BusNum),ones(1,BusNum),zeros(1,BusNum + 2 * GenNum + 2 * BrNum)];
   lc=0;
  spnum = 0;
   for l=1:ldlvnum
      b(1:2*BusNum)= datab(:,l);
    
         spnum=spnum+1;
         [plc]=nor_mskopt(A,b,c);
    
      lc = lc + plc*ldlv(l,2);
  end