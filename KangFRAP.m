function y=KangFRAP(t,re,rn,D,K)
        % Refer to Kang et al. 2008 for a discussion of this equation.
        % Partial summation from m=0 to m=10
        m=10;
        a=((-K).^(0:m)'*re^2) ./factorial(0:m)'*ones(1,length(t));
        b=re^2+((0:m)'*(8*D*t+rn^2));
        c=a./b;
        y=sum(c,1);
        
    end