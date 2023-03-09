function [lc,spnum,sp]=sp_ld17_cal(mpc,ldlv,datab,sp,spnum,totsp,CtgBrList)
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
   CtgBrNum = size(CtgBrList,1);
   for l=1:ldlvnum
      b(1:2*BusNum)= datab(:,l);
      flag=0;
      for i=spnum:-1:max(spnum-11,1) 
          ii = mod(i-1,totsp)+1;
          if (sp(ii).CtgBrNum == CtgBrNum)
              if (CtgBrNum == 0)
                  if sp(ii).invB*b > -1e-8
                     plc=sp(ii).w*b;
%                      sp(ii).num = sp(ii).num + 1;
                     flag=1;
                     break;
                  end
              elseif (sp(ii).CtgBrList == CtgBrList)
                  if sp(ii).invB*b > -1e-8
                     plc=sp(ii).w*b;
%                      sp(ii).num = sp(ii).num + 1;
                     flag=1;
                     break;
                  end
              else
                  if (sum(sp(ii).xbflag(3*BusNum+2*GenNum+CtgBrList))==CtgBrNum) && (sum(sp(ii).xbflag(3*BusNum+2*GenNum+BrNum+CtgBrList))==CtgBrNum)
                          invB = inv(A(:,sp(ii).xb));
                          w = c(sp(ii).xb)* invB;
                          sigma = c(sp(ii).xn) - w*A(:,sp(ii).xn);
                          if (sigma >= 0) & (invB(sp(ii).judge,:)*b > -1e-8)
                             sp(ii).invB = invB(sp(ii).judge,:);
                             sp(ii).w = w;
                             sp(ii).CtgBrList = CtgBrList;
                             sp(ii).CtgBrNum = size(CtgBrList,1);   
%                              sp(ii).num = sp(ii).num + 1;
                             plc=w*b;
                             flag=1;
                             break;
                          end 
                  end
              end
          else
                  if (sum(sp(ii).xbflag(3*BusNum+2*GenNum+CtgBrList))==CtgBrNum) && (sum(sp(ii).xbflag(3*BusNum+2*GenNum+BrNum+CtgBrList))==CtgBrNum)
                             invB = inv(A(:,sp(ii).xb));
                             w = c(sp(ii).xb)* invB;
                             sigma = c(sp(ii).xn) - w*A(:,sp(ii).xn);

                             if (sigma >= 0) & (invB(sp(ii).judge,:)*b > -1e-8)
                                 sp(ii).invB = invB(sp(ii).judge,:);
                                 sp(ii).w = w;
                                 sp(ii).CtgBrList = CtgBrList;
                                 sp(ii).CtgBrNum = size(CtgBrList,1);   
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
                for j=i:1:spnum-1  %max(spnum-11,1)
                           tmp = sp(mod(j-1,totsp)+1);
                           sp(mod(j-1,totsp)+1) = sp(mod(j,totsp)+1);
                           sp(mod(j,totsp)+1) = tmp;
                 end
%                 for j=spnum:-1:i+1  %max(spnum-11,1)
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
      lc = lc + plc*ldlv(l,2);
  end