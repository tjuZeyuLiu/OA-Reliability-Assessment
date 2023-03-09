%% CtgBr{1}代表0阶线路故障 
%% CtgBr{2}代表1阶线路故障 
%% CtgBr{3}代表2阶线路故障 
CtgBr{2}(1:38,1) =1;  %%CtgBr2 代表1阶线路故障
CtgBr{2}(1:38,2) =CtgList{1}(33:70,1);

CtgBr{1}(1:32,1) =1;  %%CtgBr1 代表0阶线路故障
CtgBr{1}(1:32,2) =CtgList{1}(1:32,1);
num1 = 32;
num2 = 38;
num3 = 0;
for i = 1 : size(CtgList{2},1)
    if (CtgList{2}(i,1) <=32 )&&(CtgList{2}(i,2) <=32 )
        num1 = num1  + 1;
        CtgBr{1}(num1,1) =2; 
        CtgBr{1}(num1,2) =i;
    else if (CtgList{2}(i,1) >=33 )&&(CtgList{2}(i,2) >=33 )
                num3 = num3  + 1;
                CtgBr{3}(num3,1) =2;  
                CtgBr{3}(num3,2) =i;
        else
            num2 = num2 + 1;
            CtgBr{2}(num2,1) =2;  
            CtgBr{2}(num2,2) =i;
        end
    end
end
CtgBr{1}(num1+1:num1+4960,1) =3;  %%CtgBr1 代表0阶线路故障
CtgBr{1}(num1+1:num1+4960,2) =1:4960;
num1 = num1 + 4960;
CtgBr{1}(num1+1:num1+35960,1) =4;  %%CtgBr1 代表0阶线路故障
CtgBr{1}(num1+1:num1+35960,2) =1:35960;
num1 = num1 + 35960;
CtgBr{1}(num1+1:num1+201376,1) =5;  %%CtgBr1 代表0阶线路故障
CtgBr{1}(num1+1:num1+201376,2) =1:201376;
num1 = num1 + 201376;
            