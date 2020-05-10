function J1 = costs_2016(tf, OPTIONS)

ye = zeros(10001,60);


sim('Sim_NLBM_Train_2016',[0 tf],OPTIONS,[])

J1 = J1Eval(ye, 0.0106396023101377);
end
