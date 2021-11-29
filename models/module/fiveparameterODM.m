function I = fiveparameterODM(V,Iph,Io,Rs,nVth,Rsh,toli,MaxIter)
% I = FIVEPARAMETERODM(V,IPH,IO,RS,NVTH,RSH,TOL,MAXITER) - Finds the solution to the basic 
%   One-Diode-Model equation: i = Iph-io·(exp((v+i·Rs)/nVth)-1)-(v+i·Rs)/Rsh
%   Using the Lambert-W function. To be used inside ONEDIODEMODEL / ONEDIODEMODEL2. 

    if Rs == 0 
    % Reduces to 4-parameters, explicit solution
       I = Iph-Io*(exp(V/nVth)-1)-V/Rsh;
    else
    % Solve using Lambert-W function (real,positive branch)
    %   th =kR*ki*Io*exp(ki*kR*(Iph+Io+V/Rs)), where ki = Rs/nVtg, kR = Rsh/(Rs+Rsh)
    %   I = kR*(Iph+Io-V/Rsh)-W(th)/ki

        if isinf(Rsh)
            kR = 1;
        else
            kR = Rsh/(Rs+Rsh);
        end
        ki = Rs/nVth;

        logth = log(kR*ki*Io)+ki*kR*(Iph+Io+V/Rs); % th can be big, so get ln(th) first
        W = lambertwlog(logth);
        I = kR*(Iph+Io-V/Rsh)-W/ki;

        % lambertwlog ensures precision within O(eps(W)) for W, but this doesn't imply sufficient
        % precission in I. Perform a few Newton-Raphson iterations directly over I, if necessary:
        % NOTE: use of exp(log(Io)+ ...) instead of Io·exp(...) helps reduce rounding error
        
        f = @(I) Iph-exp(log(Io)+(V+I*Rs)/nVth)+Io -(V+I*Rs)/Rsh - I;
        df = @(I) -Rs/Rsh-Rs/nVth*exp(log(Io)+(V+I*Rs)/nVth) - 1;
        if any(abs(f(I)) > toli(I))
            I = quickfzero(f,df,I,toli,MaxIter);
        end
    end
end