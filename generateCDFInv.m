function samplerand = generateCDFInv(r,mu,sigma)
erfA = erf((mu+1)/(sqrt(2)*sigma));
erfB = erf((mu-1)/(sqrt(2)*sigma));
samplerand = erfinv(-erfA-r*erfB+r*erfA)*sigma*sqrt(2)+mu;
end