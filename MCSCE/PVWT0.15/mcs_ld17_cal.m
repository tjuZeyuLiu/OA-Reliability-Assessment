function lc=mcs_ld17_cal(mpc,ll)
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
   b=[mpc.bus(:,1);mpc.bus(:,1);mpc.gen(:,2);mpc.branch(:,4);mpc.branch(:,4);];
   c=[zeros(1,BusNum),ones(1,BusNum),zeros(1,BusNum + 2 * GenNum + 2 * BrNum)];
   lc=0;
   
   lc=0;
   b(1:2*BusNum+GenNum)= ll;
   lc=nor_mskopt(A,b,c);


