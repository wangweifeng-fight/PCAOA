clear
clc
fhd=str2func('cec13_func');

func_num=1;
d = 30;

Max_iteration=1000;


tic
%
[Best_score1,Best_pos1,Convergence_curve1]=CAOA(Max_iteration,fhd,d,func_num);
toc
%
hold on
plot(Convergence_curve1,'--m','LineWidth',1.5);


hold off
legend('CAOA')



%
grid on
title('F1:Convergence curve')
xlabel('Iteration');
ylabel('Fitness function value');
% print(gcf,'-depsc','f19.eps')
