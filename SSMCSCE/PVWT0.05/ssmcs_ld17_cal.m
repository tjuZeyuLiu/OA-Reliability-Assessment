function [lc,spnum,sp]=ssmcs_ld17_cal(A,b,c,FailBrList,mpc,SPLoca,LoadS,W,datab,databgen,sp,spnum,totsp,CtgBrList,spnumdec)
  BusNum =size(mpc.bus,1);
   GenNum =size(mpc.gen,1);
   BrNum=size(mpc.branch,1);
  ldlvnum = size(SPLoca,1);
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
          CtgBrNum = size(CtgBrList,1);
    BNum = size(A,1);
    lc = 0;
    for l = 1: ldlvnum
       b(1:2*BusNum)= datab(:,LoadS(SPLoca(l,1)));
       b(2*BusNum+1:2*BusNum+GenNum)= databgen(:,LoadS(SPLoca(l,1)));

      flag=0;
      for i=spnum:-1:max(spnum-spnumdec,1) 
          ii = mod(i-1,totsp)+1;
          if (sp(ii).CtgBrNum == CtgBrNum)
              if (CtgBrNum == 0)
                  if sp(ii).judgeinvB*b > -1e-8
                     plc=sp(ii).w*b;
                     flag=1;
%                      sp(ii).num = sp(ii).num + 1;
                     break;
                  end
              elseif (sp(ii).CtgBrList == CtgBrList) %%初始矩阵B一致，这样可以不用判断检验数
                  if sp(ii).judgeinvB*b > -1e-8
                     plc=sp(ii).w*b;
                     flag=1;
%                      sp(ii).num = sp(ii).num + 1;
                     break;
                  end
              else
                  if (sum(sp(ii).xbflag(3*BusNum+2*GenNum+CtgBrList))==CtgBrNum) && (sum(sp(ii).xbflag(3*BusNum+2*GenNum+BrNum+CtgBrList))==CtgBrNum)%%要保证故障线路在基中
                          invB=fastinverse(sp(ii).invB,sp(ii).B,A(:,sp(ii).xb),BNum);
%                           invB = inv(A(:,sp(ii).xb));
                          w = c(sp(ii).xb)* invB;
                          sigma = c(sp(ii).xn) - w*A(:,sp(ii).xn);
                          if (sigma >= 0) & (invB(sp(ii).judge,:)*b > -1e-8)

                                  sp(ii).invB = invB;
                                  sp(ii).B = A(:,sp(ii).xb);
                                  sp(ii).judgeinvB = invB(sp(ii).judge,:);
                                  sp(ii).w = w;
                                  sp(ii).CtgBrList = CtgBrList;
                                  sp(ii).CtgBrNum = size(CtgBrList,1);   

                             plc=w*b;
                             flag=1;
                             break;
                          end 
                  end
              end
          else
                  if (sum(sp(ii).xbflag(3*BusNum+2*GenNum+CtgBrList))==CtgBrNum) && (sum(sp(ii).xbflag(3*BusNum+2*GenNum+BrNum+CtgBrList))==CtgBrNum)
                             invB=fastinverse(sp(ii).invB,sp(ii).B,A(:,sp(ii).xb),BNum);
                             w = c(sp(ii).xb)* invB;
                             sigma = c(sp(ii).xn) - w*A(:,sp(ii).xn);

                             if (sigma >= 0) & (invB(sp(ii).judge,:)*b > -1e-8)

                                      sp(ii).invB = invB;
                                      sp(ii).B = A(:,sp(ii).xb);
                                      sp(ii).judgeinvB = invB(sp(ii).judge,:);
                                      sp(ii).w = w;
                                      sp(ii).CtgBrList = CtgBrList;
                                      sp(ii).CtgBrNum = size(CtgBrList,1);   

                                 plc=w*b;
                                 flag=1;
                                 break;
                             end  
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
      lc = lc + plc.*W(SPLoca(l,1));
    end