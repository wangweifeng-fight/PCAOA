

function [Best_FF,Best_P,Conv_curve]=CAOA(Np,M_Iter,LB,UB,Dim,F_obj)

%Two variables to keep the positions and the fitness value of the best-obtained solution

Best_P=zeros(1,Dim);
Best_FF=inf;
Conv_curve=zeros(1,M_Iter);
lamda=10;



MOP_Max=1;
MOP_Min=0.2;
Alpha=5;
Mu=0.499;
%初始化
mu=zeros(1,Dim);
sicma=lamda*ones(1,Dim);
best=zeros(1,Dim);

for i=1:Dim
    %利用概率密度函数产生x1
    best(i)=generateCDFInv(rand,mu(i),sicma(i));
end
best=best*((UB-LB)/2)+LB;
fmin=F_obj(best);
x=zeros(1,Dim);
xnew=zeros(1,Dim);


C_Iter=1;
while C_Iter<M_Iter+1  %Main loop
    MOP=1-((C_Iter)^(1/Alpha)/(M_Iter)^(1/Alpha));   % Probability Ratio
    MOA=MOP_Min+C_Iter*((MOP_Max-MOP_Min)/M_Iter); %Accelerated function
    for i=1:Dim
        x(i)=generateCDFInv(rand,mu(i),sicma(i));
    end
    x=x*((UB-LB)/2)+LB;
    for i=1:size(x,2)
        r1=rand();
        if r1<MOA
            r2=rand();
            if r2>0.5
                xnew(i)=best(i)/(MOP+eps)*((UB-LB)*Mu+LB);
            else
                xnew(i)=best(i)*MOP*((UB-LB)*Mu+LB);
            end
        else
            r3=rand();
            if r3>0.5
                xnew(i)=best(i)-MOP*((UB-LB)*Mu+LB);
            else
                xnew(i)=best(i)+MOP*((UB-LB)*Mu+LB);
            end
        end
    end
    
    tmp=xnew;
    Flag4ub=xnew>UB;
    tmp(Flag4ub)=UB;
    Flag4lb=xnew<LB;
    tmp(Flag4lb)=LB;
    xnew=tmp;
    x=xnew;
    winner=best;
    loser=xnew;
    Fnew=F_obj(xnew);
    
    if Fnew<fmin
        winner=xnew;
        loser=best;
    end
    
    %更新mu和sigma
    munew=mu;
    mu=updateMuPV(winner,loser,munew,Np,Dim);
    sicma=updateSicmaPV(winner,loser,munew,sicma,Np,Dim);
    if Fnew<fmin
        fmin=Fnew;
        best=xnew;
    end



if fmin<Best_FF
    Best_FF=fmin;
    Best_P=best;
end

%Update the convergence curve
Conv_curve(C_Iter)=Best_FF;


C_Iter=C_Iter+1;  % incremental iteration
end
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






