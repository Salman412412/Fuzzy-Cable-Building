function J1 = J1Eval(ye, d_max)
idx_disp= 3:3:size(ye,2);
hi = [ 5.4864 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 ...
       3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 3.9624 ...
       3.9624 3.9624 3.9624 3.9624 ];

di   = [ye(:,idx_disp(1)) diff(ye(:,idx_disp)')'];

J1 = max(max(abs(di))./hi) / d_max;
end