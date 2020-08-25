%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_C_k(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSE) / MSR;

%F = (1-r0)/(MSR/MSE);
F = (MSR/MSE) * (1-r0);
df1 = n - 1;
df2 = (n-1)*(k-1);
p = 1-fcdf(F, df1, df2);

FL = (MSR/MSE) / finv(1-alpha/2, n-1, (n-1)*(k-1));
FU = (MSR/MSE) * finv(1-alpha/2, (n-1)*(k-1), n-1);

LB = 1 - 1 / FL;
UB = 1 - 1 / FU;