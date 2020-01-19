%A Program to simulate the MMSE
clc;
clear variables;
close all;
% fraction_known_variance = zeros(1000,1);
% fraction_unknown_variance = zeros(1000,1);
% for l = 1:1000
%   n = 10*l;
%   E_b = 1; N_0 = 2; m = 100; B = -1;
%   bit_error_rate = 0.5*erfc(sqrt(E_b/N_0));
%   var_actual = bit_error_rate.*(1-bit_error_rate);
%   %generate the recieved values
%   rng(592,'twister');
%   Noice = randn(n,m);
%   R = sqrt(E_b)*B + sqrt(N_0/2).*Noice;
%   if B<0
%        R = -1.*R;
%   end
%   X = zeros(n,m);
%   %Measuring the error
%   for i = 1:n
%       for j = 1:m
%           if R(i,j)<0
%               X(i,j)=1;
%           end
%       end
%   end
%   %estimating the mean
%   mean_estimate = mean(X);
%   %Calculating the interval for mean with known variance
%   delta = sqrt(var_actual)/(sqrt(n));
%   mean_upper_limit_KV = mean_estimate + delta;
%   mean_lower_limit_KV = mean_estimate - delta;
%   %Calculating the interval for mean with unknown variance
%   std_estimate = std(X);
%   deltaunk = std_estimate./(sqrt(n));
%   mean_upper_limit_unk = mean_estimate + deltaunk;
%   mean_lower_limit_unk = mean_estimate - deltaunk;
%   %Counting fraction
%   count_KV = 0;
%   countunk = 0;
%   for i = 1:m
%     if bit_error_rate<=mean_upper_limit_KV(1,i) && bit_error_rate>=mean_lower_limit_KV(1,i)
%         count_KV = count_KV + 1;
%     end
%     if bit_error_rate<=mean_upper_limit_unk(1,i) && bit_error_rate>=mean_lower_limit_unk(1,i)
%         countunk = countunk + 1;
%     end
%   end
%   fraction_known_variance(l,1) = (count_KV/m);
%   fraction_unknown_variance(l,1) = (countunk/m);
% 
% end
% figure();
% p1 = histogram(fraction_known_variance,'Normalization','probability');
% l1 = ('Known Variance');
% %hold on;
% %p2 = histogram(fraction_unknown_variance);
% %l2 = ('Unknown Variance');
% %p3 = plot (0.683*ones(1,1000));
% %l3 = ('Desired Fraction');
% %legend([p1 p2], {l1, l2});
%E_b = 1; N_0 = 2; m = 100; B = -1;
E_b = input('Enter the value of power of Signal(E_b): ');
N_0 = input('Enter the value for power of Noice(N_0) : ');
n = input('Enter the number of bits: ');
m = input('Enter the number of trials: ');
B = input('Enter the bit value(restricted to -1,1): ');
while B~=-1 && B~=1
    B = input('You did not enter the right value, Please re-enter: ');
end
bit_error_rate = 0.5*erfc(sqrt(E_b/N_0));
var_actual = bit_error_rate.*(1-bit_error_rate);
%generate the recieved values
rng(592,'twister');
Noice = randn(n,m);
R = sqrt(E_b)*B + sqrt(N_0/2).*Noice;
if B<0
    R = -1.*R;
end
X = zeros(n,m);
%Measuring the error
for i = 1:n
    for j = 1:m
        if R(i,j)<0
            X(i,j)=1;
        end
    end
end
%estimating the mean
mean_estimate = mean(X);
%Calculating the interval for mean with known variance
delta = sqrt(var_actual)/(sqrt(n));
mean_upper_limit_KV = mean_estimate + delta;
mean_lower_limit_KV = mean_estimate - delta;
%Calculating the interval for mean with unknown variance
std_estimate = std(X);
deltaunk = std_estimate./(sqrt(n));
mean_upper_limit_unk = mean_estimate + deltaunk;
mean_lower_limit_unk = mean_estimate - deltaunk;
%Counting fraction
count_KV = 0;
countunk = 0;
for i = 1:m
    if bit_error_rate<=mean_upper_limit_KV(1,i) && bit_error_rate>=mean_lower_limit_KV(1,i)
        count_KV = count_KV + 1;
    end
    if bit_error_rate<=mean_upper_limit_unk(1,i) && bit_error_rate>=mean_lower_limit_unk(1,i)
        countunk = countunk + 1;
    end
end
fraction_known_variance = (count_KV/m);
fraction_unknown_variance = (countunk/m);
% %Plotting all the outcomes
figure();
p1=plot(1:100, mean_estimate(1:100));
l1='Estimated BER';
hold on
p2=plot(1:100, mean_lower_limit_KV(1:100));
l2='Lower Level of Confidence Interval(68.3%)(Known variance)';
p3=plot(1:100, mean_upper_limit_KV(1:100));
l3='Upper Level of Confidence Interva(68.3%)(Known variance)';
p4=plot(1:100, bit_error_rate*ones(1,100));
l4='True bit error rate (BER) or probability of error';
hold off
legend([p1 p2 p3 p4], {l1, l2, l3, l4});
xlabel('Trials');
title('Comparision of Estimated Mean vs True Mean(True Variance)');
figure();
p5=plot(1:100, mean_estimate(1:100));
l5='Estimated Value of BER';
hold on
p6=plot(1:100, mean_lower_limit_unk(1:100));
l6='Lower Level of Confidence Interval(68.3%)(estimated variance)';
p7=plot(1:100, mean_upper_limit_unk(1:100));
l7='Upper Level of Confidence Interval(68.3%)(estimated variance)';
p8=plot(1:100, bit_error_rate*ones(1,100));
l8='True bit error rate (BER)';
hold off
legend([p5 p6 p7 p8], {l5, l6, l7, l8});
xlabel('Trials');
title('Comparision of Estimated Mean vs True Mean(Estimated Variance)');
