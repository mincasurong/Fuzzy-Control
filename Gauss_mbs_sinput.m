function [mbsfn, mbsfn_sum] = Gauss_mbs_sinput(x,c,N,s)

m=2.5;
mbsfn = zeros(N,length(c));
mbsfn_sum=zeros(N,1);
cvar=(c(2)-c(1))/m;
% eqn = @(syms_s) exp(-0.5*abs(cvar/syms_s)^m) == 0.5;
% s=round(double(solve(eqn)),4);

for k=1:N,
    for j=1:length(c), mbsfn(k,j) = round(exp(-0.5*abs((x(k)-c(j))/s)^m),3); end
    if x(k) < c(1), mbsfn(k,1) = 1; end
    if x(k) > c(length(c)), mbsfn(k,length(c)) = 1; end
    mbsfn_sum(k) = sum(mbsfn(k,:));
end