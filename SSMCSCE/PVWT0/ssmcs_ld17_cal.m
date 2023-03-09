function [lc,spnum,sp]=ssmcs_ld17_cal(mpc,SPLoca,LoadS,W,datab,sp,spnum,totsp,CtgBrList)
  BusNum =size(mpc.bus,1);
   GenNum =size(mpc.gen,1);
   BrNum=size(mpc.branch,1);
  ldlvnum = size(SPLoca,1);
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
   b=[zeros(2*BusNum,1);mpc.gen(:,2);mpc.branch(:,4);mpc.branch(:,4);];
   c=[zeros(1,BusNum),ones(1,BusNum),zeros(1,BusNum + 2 * GenNum + 2 * BrNum)];
          CtgBrNum = size(CtgBrList,1);
    BNum = size(A,1);
    lc = 0;
    for j = 1: ldlvnum
      b(1:2*BusNum)= datab(:,LoadS(SPLoca(j,1)));

      flag=0;
      for i=spnum:-1:max(spnum-20,1) 
          ii = mod(i-1,totsp)+1;
          if (sp(ii).CtgBrNum == CtgBrNum)
              if (CtgBrNum == 0)
                  if sp(ii).judgeinvB*b > -1e-8
                     plc=sp(ii).w*b;
                     flag=1;
%                      sp(ii).num = sp(ii).num + 1;
                     break;
                  end
              elseif (sp(ii).CtgBrList == CtgBrList)
                  if sp(ii).judgeinvB*b > -1e-8
                     plc=sp(ii).w*b;
                     flag=1;
%                      sp(ii).num = sp(ii).num + 1;
                     break;
                  end
              else
                  if (sum(sp(ii).xbflag(3*BusNum+2*GenNum+CtgBrList))==CtgBrNum) && (sum(sp(ii).xbflag(3*BusNum+2*GenNum+BrNum+CtgBrList))==CtgBrNum)
                          invB=fastinverse(sp(ii).invB,sp(ii).B,A(:,sp(ii).xb),BNum);
%                           invB = inv(A(:,sp(ii).xb));
                          w = c(sp(ii).xb)* invB;
                          sigma = c(sp(ii).xn) - w*A(:,sp(ii).xn);
                          if (sigma >= 0) & (invB(sp(ii).judge,:)*b > -1e-8)
%                              sp(ii).invB = invB(sp(ii).judge,:);
%                              sp(ii).w = w;
%                              sp(ii).CtgBrList = CtgBrList;
%                              sp(ii).CtgBrNum = size(CtgBrList,1);   
%                              sp(ii).num = sp(ii).num + 1;
                             plc=w*b;
                             flag=1;
                             break;
                          end 
                  end
              end
          else
                  if (sum(sp(ii).xbflag(3*BusNum+2*GenNum+CtgBrList))==CtgBrNum) && (sum(sp(ii).xbflag(3*BusNum+2*GenNum+BrNum+CtgBrList))==CtgBrNum)
                             invB=fastinverse(sp(ii).invB,sp(ii).B,A(:,sp(ii).xb),BNum);
%                                invB = inv(A(:,sp(ii).xb));
                             w = c(sp(ii).xb)* invB;
                             sigma = c(sp(ii).xn) - w*A(:,sp(ii).xn);

                             if (sigma >= 0) & (invB(sp(ii).judge,:)*b > -1e-8)
%                                  sp(ii).invB = invB(sp(ii).judge,:);
%                                  sp(ii).w = w;
%                                  sp(ii).CtgBrList = CtgBrList;
%                                  sp(ii).CtgBrNum = size(CtgBrList,1); 
%                                  sp(ii).num = sp(ii).num + 1;
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
%                if sp(ii).num > sp(mod(i+1-1,totsp)+1).num
%                    tmp = sp(ii);
%                    sp(ii) = sp(mod(i+1-1,totsp)+1);
%                    sp(mod(i+1-1,totsp)+1) = tmp;
%                end
                for kkk=i:1:spnum-1  %max(spnum-11,1)
                           tmp = sp(mod(kkk-1,totsp)+1);
                           sp(mod(kkk-1,totsp)+1) = sp(mod(kkk,totsp)+1);
                           sp(mod(kkk,totsp)+1) = tmp;
                 end
%                tmp = sp(ii);
%                sp(ii) = sp(mod(spnum-1,totsp)+1);
%                sp(mod(spnum-1,totsp)+1) = tmp;
               
%                                for j=i:1:spnum-1  %max(spnum-11,1)
%                            tmp = sp(mod(j-1,totsp)+1);
%                            sp(mod(j-1,totsp)+1) = sp(mod(j,totsp)+1);
%                            sp(mod(j,totsp)+1) = tmp;
%                  end
%                end
%                  for j=spnum:-1:i+1  %max(spnum-11,1)
%                        if sp(mod(j-1,totsp)+1).num <= sp(ii).num
%                            
%                            for k = i : 1 : j-1
%                                    tmp = sp(mod(k-1,totsp)+1);
%                                    sp(mod(k-1,totsp)+1) = sp(mod(k,totsp)+1);
%                                    sp(mod(k,totsp)+1) = tmp;
%                            end
%                            break;
%                        end
%                  end
          end               
      end
      lc = lc + plc.*W(SPLoca(j,1));
    end

