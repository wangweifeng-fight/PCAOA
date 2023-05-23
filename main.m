%_______________________________________________________________________________________
%  The Arithmetic Optimization Algorithm (AOA) source codes demo version 1.0                  
%_______________________________________________________________________________________

clear all 
clc


% Solution_no=20; %Number of search solutions
F_name='F1';    %Name of the test function F1-f23
M_Iter=500;    %Maximum number of iterations
 
[LB,UB,Dim,F_obj]=Get_F(F_name); %Give details of the underlying benchmark function
Solution_no=2*Dim;
[Best_FF,Best_P,Conv_curve]=AOA(Solution_no,M_Iter,LB,UB,Dim,F_obj); % Call the AOA 

 

figure('Position',[454   445   694   297]);
subplot(1,2,1);
func_plot(F_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([F_name,'( x_1 , x_2 )'])


subplot(1,2,2);
semilogy(Conv_curve,'Color','r','LineWidth',2)
title('Convergence curve')
xlabel('Iteration#');
ylabel('Best fitness function');
axis tight
legend('AOA')



display(['The best-obtained solution by Math Optimizer is : ', num2str(Best_P)]);
display(['The best optimal value of the objective funciton found by Math Optimizer is : ', num2str(Best_FF)]);

        



