function [num,sum_num,mu] = fuzzyrule_specific(a,b,u_mbs0)

num = zeros(15,1);
mu_num = zeros(15,1);

num(1) = a(1) * b(1); mu_num(1) = num(1) * u_mbs0(1); % NH
num(2) = a(2) * b(1); mu_num(2) = num(2) * u_mbs0(2);  % NL
num(3) = a(3) * b(1); mu_num(3) = num(3) * u_mbs0(3);  % ZN
num(4) = a(4) * b(1); mu_num(4) = num(4) * u_mbs0(2);  % NL
num(5) = a(5) * b(1); mu_num(5) = num(5) * u_mbs0(1);  % NH

num(6) = a(1) * b(2); mu_num(6) = num(6) * u_mbs0(2);  % NL
num(7) = a(2) * b(2); mu_num(7) = num(7) * u_mbs0(3);  % ZN
num(8) = a(3) * b(2); mu_num(8) = num(8) * 0;          % ZE
num(9) = a(4) * b(2); mu_num(9) = num(9) * u_mbs0(4);  % ZP
num(10) = a(5) * b(2);mu_num(10) = num(10) * u_mbs0(5);   % PL

num(11) = a(1) * b(3); mu_num(11) = num(11) * u_mbs0(6);  % PH
num(12) = a(2) * b(3); mu_num(12) = num(12) * u_mbs0(5);  % PL
num(13) = a(3) * b(3); mu_num(13) = num(13) * u_mbs0(4);  % ZP
num(14) = a(4) * b(3); mu_num(14) = num(14) * u_mbs0(5);  % PL
num(15) = a(5) * b(3); mu_num(15) = num(15) * u_mbs0(6);  % PH

sum_num=sum(num); 
mu = sum(mu_num)/sum_num;