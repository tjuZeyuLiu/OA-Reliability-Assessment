function [lc,spnum,sp]=ssmcs_ld17_cal_init(mpc,SPLoca,LoadS,W,datab,totsp,CtgBrList)
   BusNum =size(mpc.bus,1);
   GenNum =size(mpc.gen,1);
   BrNum=size(mpc.branch,1);
   ldlvnum = size(SPLoca,1);
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
       c=[zeros(1,BusNum),ones(1,BusNum),zeros(1,BusNum + 2 * GenNum + 2 * BrNum)];

       sp(1:totsp) = struct('A',[],'B',[],'invB',[],'judgeinvB',[],'w',[],'xb',[],'xn',[],'CtgBrList',[],'CtgBrNum',[],'judge',[],'xbflag',[],'num',[]);
       spnum=0;
       CtgBrNum = size(CtgBrList,1);
  lc = 0;
    for j = 1: ldlvnum

       b(1:2*BusNum)= datab(:,LoadS(SPLoca(j,1)));

      flag=0;
      for i=spnum:-1:max(spnum-4,1)
          if ((sp(mod(i-1,totsp)+1).CtgBrNum == 0)&& (CtgBrNum == 0))|(sp(mod(i-1,totsp)+1).CtgBrList == CtgBrList)
          if sp(mod(i-1,totsp)+1).invB*b > -1e-8
            plc = sp(mod(i-1,totsp)+1).w*b;
%              sp(mod(i-1,totsp)+1).num = sp(mod(i-1,totsp)+1).num + 1;
             flag=1;
             break;
          end
          end
      end
      if flag==0
         spnum=spnum+1;
         [plc,sp(mod(spnum-1,totsp)+1)]=sp_newmskopt(A,b,c,BusNum,GenNum,CtgBrList);
       else
          if i~=spnum

                for kkk=i:1:spnum-1  %max(spnum-11,1)
                           tmp = sp(mod(kkk-1,totsp)+1);
                           sp(mod(kkk-1,totsp)+1) = sp(mod(kkk,totsp)+1);
                           sp(mod(kkk,totsp)+1) = tmp;
                 end

          end
      end
      lc = lc + plc.*W(SPLoca(j,1));
    end
