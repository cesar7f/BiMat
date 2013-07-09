ssm = mod.ss;
ss = zeros(34,4);
for i = 1:34
    ss(i,ssm(i)) = 1;
end