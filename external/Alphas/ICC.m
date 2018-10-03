function [r, LB, UB, F, df1, df2, p] = ICC(M, type, alpha, r0)
% Intraclass correlation
%   [r, LB, UB, F, df1, df2, p] = ICC(M, type, alpha, r0)
%
%   M is matrix of observations. Each row is an object of measurement and
%   each column is a judge or measurement.
%
%   'type' is a string that can be one of the six possible codes for the desired
%   type of ICC:
%       '1-1': The degree of absolute agreement among measurements made on
%         randomly seleted objects. It estimates the correlation of any two
%         measurements.
%       '1-k': The degree of absolute agreement of measurements that are
%         averages of k independent measurements on randomly selected
%         objects.
%       'C-1': case 2: The degree of consistency among measurements. Also known
%         as norm-referenced reliability and as Winer's adjustment for
%         anchor points. case 3: The degree of consistency among measurements maded under
%         the fixed levels of the column factor. This ICC estimates the
%         corrlation of any two measurements, but when interaction is
%         present, it underestimates reliability.
%       'C-k': case 2: The degree of consistency for measurements that are
%         averages of k independent measurements on randomly selected
%         onbjectgs. Known as Cronbach's alpha in psychometrics. case 3:  
%         The degree of consistency for averages of k independent
%         measures made under the fixed levels of column factor.
%       'A-1': case 2: The degree of absolute agreement among measurements. Also
%         known as criterion-referenced reliability. case 3: The absolute 
%         agreement of measurements made under the fixed levels of the column factor.
%       'A-k': case 2: The degree of absolute agreement for measurements that are
%         averages of k independent measurements on randomly selected objects.
%         case 3: he degree of absolute agreement for measurements that are
%         based on k independent measurements maded under the fixed levels
%         of the column factor.
%
%       ICC is the estimated intraclass correlation. LB and UB are upper
%       and lower bounds of the ICC with alpha level of significance. 
%
%       In addition to estimation of ICC, a hypothesis test is performed
%       with the null hypothesis that ICC = r0. The F value, degrees of
%       freedom and the corresponding p-value of the this test are
%       reported.
%
%       (c) Arash Salarian, 2008
%
%       Reference: McGraw, K. O., Wong, S. P., "Forming Inferences About
%       Some Intraclass Correlation Coefficients", Psychological Methods,
%       Vol. 1, No. 1, pp. 30-46, 1996
%

if nargin < 3
    alpha = .05;
end

if nargin < 4
    r0 = 0;
end

[n, k] = size(M);

SStotal = var(M(:)) *(n*k - 1);
MSR = var(mean(M, 2)) * k;
MSW = sum(var(M,0, 2)) / n;
MSC = var(mean(M, 1)) * n;
MSE = (SStotal - MSR *(n - 1) - MSC * (k -1))/ ((n - 1) * (k - 1));

switch type
    case '1-1'
        [r, LB, UB, F, df1, df2, p] = ICC_case_1_1(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case '1-k'
        [r, LB, UB, F, df1, df2, p] = ICC_case_1_k(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case 'C-1'
        [r, LB, UB, F, df1, df2, p] = ICC_case_C_1(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case 'C-k'
        [r, LB, UB, F, df1, df2, p] = ICC_case_C_k(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case 'A-1'
        [r, LB, UB, F, df1, df2, p] = ICC_case_A_1(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case 'A-k'
        [r, LB, UB, F, df1, df2, p] = ICC_case_A_k(MSR, MSE, MSC, MSW, alpha, r0, n, k);
end