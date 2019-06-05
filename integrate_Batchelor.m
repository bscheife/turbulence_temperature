
function [I]=integrate_Batchelor(k1,k2,kB,X_T,Dt,q)
%Calculate the integral of the Batchelor spectrum between k1 and k2, given
%kB, Chi, Dt, and q.

%Batchelor spectral shape
funcn=@(alpha)(alpha.*(exp(-alpha.^2/2)-alpha.*sqrt(pi/2).*erfc(alpha/sqrt(2))));
%Integration limits, in alpha space
alpha1=sqrt(2*q)*k1/kB;
alpha2=sqrt(2*q)*k2/kB;
%Integrate the above function from alpha1 to alpha2
I = X_T/(2*Dt) * quad(funcn,alpha1,alpha2);

end



%Verified B.Scheifele 2017-06: 
%int(SB)dk from k1 to k2 is equivalent to X/(2*Dt)*int(funcn)d{alpha} from
%alpha1 to alpha2, where funcn is as defined above and SB is the Batchelor
%spectrum given by:
    % alpha=sqrt(2*q)*(k/kB);
    % SB = sqrt(q/2) * X_T ./ (kB.*Dt) .* ...
    %     alpha.*(exp(-alpha.^2/2)-sqrt(pi/2)*alpha.*erfc(alpha/sqrt(2)));
%as coded in Batch_spec.m
