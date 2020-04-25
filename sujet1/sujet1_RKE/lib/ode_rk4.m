function  [T,X] = ode_rk4(f,t0tf,x0,N)
%-------------------------------------------------------------------------------------------
%
%    RK4 method 
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
%        [T,X] = ode_rk4(f,t0tf,y0,N)
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
     
      h = (t0tf(2)-t0tf(1))/N;
     T(1) = 0;
     X(1,:) = x0 ; 
     
     for i = 2:N+1  
         
         T(i) = (i-1)*h;
         k1 = f(T(i-1),X(i-1,:))';
         k2 = f(T(i-1) + (1/2)*h, X(i-1,:) + h*(1/2)*k1)';  
         k3 = f(T(i-1) + (1/2)*h, X(i-1,:) + h*(0*k1 + (1/2)*k2))';
         k4 = f(T(i-1) + h,X(i-1,:) + h*(0*k1 + 0*k2 + k3))';
         X(i,:) = X(i-1,:) + h*((1/6)*k1 +(2/6)*k2 + (2/6)*k3 + (1/6)*k4 );         
         
         
     end 
     
end


