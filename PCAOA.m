
function [Best_FF,Best_P,Conv_curve]=PCAOA(Np,M_Iter,LB,UB,Dim,F_obj)

%Two variables to keep the positions and the fitness value of the best-obtained solution

Best_P=zeros(1,Dim);
Best_FF=inf;
Conv_curve=zeros(1,M_Iter);
lamda=10;
groups=4;

Group=repmat(struct('mu',{},'sigma',{}),1,groups);

%初始化
for g=1:groups
    Group(g).mu=zeros(1,Dim);
    Group(g).sicma=lamda*ones(1,Dim);
    Group(g).best=zeros(1,Dim);
    for i=1:Dim
        %利用概率密度函数产生x1
        Group(g).best(i)=generateCDFInv(rand,Group(g).mu(i),Group(g).sicma(i));
    end
    Group(g).best=Group(g).best*((UB-LB)/2)+(UB+LB)/2;
    Group(g).fmin=F_obj(Group(g).best);
    Group(g).x=zeros(1,Dim);
    Group(g).xnew=zeros(1,Dim);
end
C_Iter=1;
while C_Iter<M_Iter+1  %Main loop
    for g=1:groups
        Group(g).xnew=aoa(Group(g).best,Group(g).x,C_Iter,M_Iter,Group,g,Dim,UB,LB);
        winner=Group(g).best;
        loser=Group(g).xnew;
        Fnew=F_obj(Group(g).xnew);
        
        if Fnew<Group(g).fmin
            winner=Group(g).xnew;
            loser=Group(g).best;
        end
        
        %更新mu和sigma
        munew=Group(g).mu;
        Group(g).mu=updateMuPV(winner,loser,munew,Np,Dim);
        Group(g).sicma=updateSicmaPV(winner,loser,munew,Group(g).sicma,Np,Dim);
        if Fnew<Group(g).fmin
            Group(g).fmin=Fnew;
            Group(g).best=Group(g).xnew;
        end
    end
    
    for g=1:groups
        if Group(g).fmin<Best_FF
            Best_FF=Group(g).fmin;
            Best_P=Group(g).best;
        end
    end
    
    for g=1:groups
        g_r1=randperm(groups,1);
        g_r2=randperm(groups,1);
        g_r3=randperm(groups,1);
        if rem(M_Iter,10)==0
            %当前组与随机组
            if Group(g).fmin > Group(g_r1).fmin
                Group(g).fmin  = Group(g_r1).fmin;
                Group(g).best = Group(g_r1).best;
            else
                average_curr=mean(Group(g).xnew);
                average_rd=mean(Group(g_r1).xnew);
                sort(Group(g).xnew,'descend');
                sort(Group(g_r1).xnew,'descend');
                if mod(Dim,2)==0
                    sorted_x_cur=Group(g).xnew(((Dim/2)+1):end);
                    sorted_x_rd=Group(g_r1).xnew(((Dim/2)+1):end);
                else
                    sorted_x_cur=Group(g).xnew((floor(Dim/2)+1):end);
                    sorted_x_rd=Group(g_r1).xnew((ceil(Dim/2)+1):end);
                end
                if average_curr>average_rd
                    new_xrd=[sorted_x_cur,sorted_x_rd];
                    fit_rd = F_obj(new_xrd);
                    if fit_rd < Group(g).fmin
                        Group(g).fmin = fit_rd;
                        Group(g).best=new_xrd;
                    end
                end
            end
            if rem(M_Iter,20)==0
                if Group(g_r2).fmin > Group(g_r3).fmin
                    Group(g_r2).fmin  = Group(g_r3).fmin;
                    Group(g_r2).best = Group(g_r3).best;
                else
                    average_curr=mean(Group(g_r2).xnew);
                    average_rd=mean(Group(g_r3).xnew);
                    sort(Group(g_r2).xnew,'descend');
                    sort(Group(g_r3).xnew,'descend');
                    if mod(Dim,2)==0
                        sorted_x_cur=Group(g_r2).xnew(((Dim/2)+1):end);
                        sorted_x_rd=Group(g_r3).xnew(((Dim/2)+1):end);
                    else
                        sorted_x_cur=Group(g_r2).xnew((floor(Dim/2)+1):end);
                        sorted_x_rd=Group(g_r3).xnew((ceil(Dim/2)+1):end);
                    end
                    if average_curr>average_rd
                        new_xrd=[sorted_x_cur,sorted_x_rd];
                        fit_rd = F_obj(new_xrd);
                        if fit_rd < Group(g_r2).fmin
                            Group(g_r2).fmin = fit_rd;
                            Group(g_r2).best=new_xrd;
                        end
                    else
                        new_xrd=[sorted_x_rd,sorted_x_cur];
                        fit_rd = F_obj(new_xrd);
                        if fit_rd < Group(g_r3).fmin
                            Group(g_r3).fmin = fit_rd;
                            Group(g_r3).best=new_xrd;
                        end
                    end                  
                end
            end
        end
        %         if rand<0.5
        %             %循环当前组和随机组、最佳组进行比较
        %             if rand<0.5
        %
        %             else
        %                 %当前组与最佳组
        %                 if Group(g).fmin > Best_FF
        %                     Group(g).fmin = Best_FF;
        %                     Group(g).best = Best_P;
        %                 else
        %                     average_curr=mean(Group(g).xnew);
        %                     average_best=mean(Best_P);
        %                     sort(Group(g).xnew,'descend');
        %                     sort(Best_P,'descend');
        %                     if mod(Dim,2)==0
        %                         sorted_x_cur=Group(g).x(((Dim/2)+1):end);
        %                         sorted_x_best=Best_P(((Dim/2)+1):end);
        %                     else
        %                         sorted_x_cur=Group(g).x((floor(Dim/2)+1):end);
        %                         sorted_x_best=Best_P((ceil(Dim/2)+1):end);
        %                     end
        %                     if average_curr>average_best
        %                         new_xrd=[sorted_x_cur,sorted_x_best];
        %                         fit_rd = F_obj(new_xrd);
        %                         if fit_rd <Group(g).fmin
        %                             Best_FF = fit_rd;
        %                             Best_P = new_xrd;
        %                         end
        %                     end
        %                 end
        %             end
        %         else
        %             %循环随机组和最佳组进行比较
        %             if Group(g_r2).fmin>Best_FF
        %                 Group(g_r2).fmin=Best_FF;
        %                 Group(g_r2).best=Best_P;
        %             else
        %                 average_rd1=mean(Group(g_r2).x);
        %                 average_best=mean(Best_P);
        %                 sort(Group(g_r2).x);
        %                 sort(Best_P,'descend');
        %                 if mod(Dim,2)==0
        %                     sort_x_rd1=Group(g_r2).x(((Dim/2)+1):end);
        %                     sort_x_best=Best_P(((Dim/2)+1):end);
        %                 else
        %                     sort_x_rd1=Group(g_r2).x((floor(Dim/2)+1):end);
        %                     sort_x_best=Best_P((ceil(Dim/2)+1):end);
        %                 end
        %                 if average_rd1>average_best
        %                     new_xrd=[sort_x_rd1,sort_x_best];
        %                     fit_rd = F_obj(new_xrd);
        %                     if fit_rd <Group(g_r2).fmin
        %                         Best_FF = fit_rd;
        %                         Best_P = new_xrd;
        %                     end
        %                 end
        %             end
        %         end
        [tmp,index ] = sort(Group(g).fmin,'descend');
        [tmp1,index1] = sort(Group(g).fmin);
        
        pos1 = index(1);%最差的个体
        pos2 = index1(1);%最优的个体
        Group(pos1).best=Group(pos2).best;
        Group(pos1).fmin=Group(pos2).fmin;
        
    end
    
    %Update the convergence curve
    Conv_curve(C_Iter)=Best_FF;
    
    
    C_Iter=C_Iter+1;  % incremental iteration
    
end
end

function x_new=aoa(best,X,C_Iter,M_Iter,Group,g,Dim,UB,LB)

MOP_Max=1;
MOP_Min=0.2;
Alpha=5;
Mu=0.499;
MOP=1-((C_Iter)^(1/Alpha)/(M_Iter)^(1/Alpha));   % Probability Ratio
MOA=MOP_Min+C_Iter*((MOP_Max-MOP_Min)/M_Iter); %Accelerated function
for i=1:Dim
    X(i)=generateCDFInv(rand,Group(g).mu(i),Group(g).sicma(i));
end
X=X*((UB-LB)/2)+(UB+LB)/2;
Xmin=min(abs(X));
Best=repmat(Xmin,1,Dim);
r1=rand();
if r1<MOA
    r2=rand();
    if r2>0.5
        x_new=Best/(MOP+eps)*((UB-LB)*Mu+LB);
    else
        x_new=Best*MOP*((UB-LB)*Mu+LB);
    end
else
    r3=rand();
    if r3>0.5
        x_new=Best-MOP*((UB-LB)*Mu+LB);
    else
        x_new=Best+MOP*((UB-LB)*Mu+LB);
    end
end
Flag_UB=x_new>UB; % check if they exceed (up) the boundaries
Flag_LB=x_new<LB; % check if they exceed (down) the boundaries
x_new=(x_new.*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;

end
function samplerand = generateCDFInv(r,mu,sigma)
erfA = erf((mu+1)/(sqrt(2)*sigma));
erfB = erf((mu-1)/(sqrt(2)*sigma));
samplerand = erfinv(-erfA-r*erfB+r*erfA)*sigma*sqrt(2)+mu;
end
function mmu =updateMuPV(win,lost,mu,Np,dim) %update meanVector belong [-1,1];
mmu=mu;
for k=1:dim
    dm =(1/Np)*(win(k)-lost(k));
    %     mmu(k)=mu(k)+dm;
    if (abs(mu(k)+dm)<=1)
        mmu(k)=mu(k)+dm;
    end
end
end
function ssicma =updateSicmaPV(win,los,mu,sicma,Np,dim)
sicma2=sicma.^2;
ssicma2=sicma2;
mmu=mu;
for k=1:dim
    dm   =(1/Np)*(win(k)-los(k));
    mmu(k)=mu(k)+dm;
    A=mu(k).^2;
    B=mmu(k).^2;
    C=(1/Np)*(win(k).^2-los(k).^2);
    dsicma2=A-B+C;
    if abs(dsicma2)<sicma2(k)
        ssicma2(k)=sicma2(k)+dsicma2;
    else
        ssicma2(k)=sicma2(k);
    end
end
ssicma=sqrt(ssicma2);
end



