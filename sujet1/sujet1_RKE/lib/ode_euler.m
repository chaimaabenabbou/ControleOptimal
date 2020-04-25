function  [T,X] = ode_euler(f,t0tf,x0,N)
%-------------------------------------------------------------------------------------------
%
%    Explicit Euler method 
%
%    Description
%
%        Numerical integration of the Cauchy's problem
%        \dot{x}(t) = f(t,x(t))
%        x(t_0) = x_0
%
%-------------------------------------------------------------------------------------------
%
%    Usage
%
%        [T,X] = ode_euler(f,t0tf,y0,N)
%
%    Inputs
%        f    - function     : second member of the ode whith the interface 
%                              xpoint = f(t, x)
%                                  t    - real     : time,
%                                  x = vector of R^n with the same dimension of x0
%        t0tf - real(2)      : intial and final time  [t0,tf]
%        x0   - real(n)      : initial point
%        N    - integer      : number of steps (>1)
%
%    Outputs
%        T    - real(N+1,1)  : vector of times
%        X    - real(N+1,n)  : Matrix of solution
%        The line i of [T Y] contains ti and x_i
%
%-------------------------------------------------------------------------------------------
     n = length(x0);
% A remplacer     
     T = zeros(N+1,1);
     X = zeros(N+1,n);
     
     pas = (t0tf(2)-t0tf(1))/N;
      T(1) = 0;
     X(1,:) = x0 ; 
    
     
     for i = 2:N+1         
         X(i,:) = X(i-1,:) + pas*f(T(i-1),X(i-1,:))';
         T(i) = (i-1)*pas;
     end 
     
     
     
end



