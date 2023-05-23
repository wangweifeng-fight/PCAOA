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