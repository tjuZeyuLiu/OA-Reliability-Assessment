a1 = 0;
a2 = 1;
t = 0.1;
for i = 1 : 10
    level(:,i) = ldlv(:,1)*(a1+(i-1)*t)+ldlv(:,2)*(a2-(i-1)*t);
end
level(:,11) =1;


save('ld11_lv8760.mat','level');