function [lc,spnum,sp]=ssmcs_ld17_cal_init(A,b,c,FailBrList,mpc,SPLoca,LoadS,W,datab,databgen,totsp,CtgBrList,spnumdec)
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
     for  i = 1 : size(FailBrList,1)
        A(mpc.branch(FailBrList(i),1),mpc.branch(FailBrList(i),2)) = A(mpc.branch(FailBrList(i),1),mpc.branch(FailBrList(i),2)) - mpc.branch(FailBrList(i),3);
        A(mpc.branch(FailBrList(i),2),mpc.branch(FailBrList(i),1)) = A(mpc.branch(FailBrList(i),2),mpc.branch(FailBrList(i),1)) - mpc.branch(FailBrList(i),3);
        A(mpc.branch(FailBrList(i),1),mpc.branch(FailBrList(i),1)) = A(mpc.branch(FailBrList(i),1),mpc.branch(FailBrList(i),1)) + mpc.branch(FailBrList(i),3);
        A(mpc.branch(FailBrList(i),2),mpc.branch(FailBrList(i),2)) = A(mpc.branch(FailBrList(i),2),mpc.branch(FailBrList(i),2)) + mpc.branch(FailBrList(i),3);
        A(GenNum+2*BusNum+FailBrList(i),mpc.branch(FailBrList(i),1)) = A(GenNum+2*BusNum+FailBrList(i),mpc.branch(FailBrList(i),1)) - mpc.branch(FailBrList(i),3);
        A(GenNum+2*BusNum+FailBrList(i),mpc.branch(FailBrList(i),2)) = A(GenNum+2*BusNum+FailBrList(i),mpc.branch(FailBrList(i),2)) + mpc.branch(FailBrList(i),3);
        A(GenNum+2*BusNum+BrNum+FailBrList(i),mpc.branch(FailBrList(i),1)) = A(GenNum+2*BusNum+BrNum+FailBrList(i),mpc.branch(FailBrList(i),1)) + mpc.branch(FailBrList(i),3);
        A(GenNum+2*BusNum+BrNum+FailBrList(i),mpc.branch(FailBrList(i),2)) = A(GenNum+2*BusNum+BrNum+FailBrList(i),mpc.branch(FailBrList(i),2)) - mpc.branch(FailBrList(i),3);
        b(2*BusNum+GenNum+FailBrList(i)) = 0;
        b(2*BusNum+GenNum+BrNum+FailBrList(i)) = 0;
     end

  sp(1:totsp) = struct('A',[],'B',[],'invB',[],'judgeinvB',[],'w',[],'xb',[],'xn',[],'CtgBrList',[],'CtgBrNum',[],'judge',[],'xbflag',[],'num',[]);
       spnum=0;
       CtgBrNum = size(CtgBrList,1);
  lc = 0;
    for j = 1: ldlvnum

       b(1:2*BusNum)= datab(:,LoadS(SPLoca(j,1)));
       b(2*BusNum+1:2*BusNum+GenNum)= databgen(:,LoadS(SPLoca(j,1)));

      flag=0;
      for i=spnum:-1:max(spnum-spnumdec,1)
          if ((sp(mod(i-1,totsp)+1).CtgBrNum == 0)&& (CtgBrNum == 0))|(sp(mod(i-1,totsp)+1).CtgBrList == CtgBrList)
          if sp(mod(i-1,totsp)+1).judgeinvB*b > -1e-8
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

                for k=i:1:spnum-1  %max(spnum-11,1)
                           tmp = sp(mod(k-1,totsp)+1);
                           sp(mod(k-1,totsp)+1) = sp(mod(k,totsp)+1);
                           sp(mod(k,totsp)+1) = tmp;
                 end

          end    
      end
      lc = lc + plc.*W(SPLoca(j,1));
    end
