function [num,sum_num,mu] = fuzzyrule(a,b,u_mbs0)

END = 5;
num = zeros(1,END);
mu_num = zeros(1,END);

num(1,1) = a(1) * b(1);  % NH
num(1,2) = a(2) * b(1);  % NL
num(1,3) = a(3) * b(2);  % ZN
num(1,4) = a(4) * b(2);  % ZP
num(1,5) = a(5) * b(3);  % PH



sum_num=sum(num); 
for k=1:END, mu_num(k)=num(k)*u_mbs0(k); end
mu = sum(mu_num)/sum_num;