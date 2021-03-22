%These are the parameters for exercise 3.
%They should be substitued into execise 2
%https://www.blogs.uni-mainz.de/fb09-geosciences/files/2018/12/Finite_Differerence_Timpl_expl.pdf

clear

% Physical parameters
L       =   100;         % Length of modeled domain     [m]
Tmagma  =   1200;        % Temperature of magma         [C]
Trock   =   300;         % Temperature of country rock   [C]
Tleft   =   200;         % Temperature of left boundary [C]
Tright  =   100;         % Temperature of right boundary [C]
kappa   =   1e-6;        % Thermal diffusivity of rock  [m2/s]
W       =   5;           % Width of dike                [m]
day     =   3600*24;     % # seconds per day
dt      =   1*day;       % Timestep                     [s]


% Numerical parameters
nx      =   201;         % Number of gridpoints in x-direction
nt      =   500;         % NUmber of timesteps to compute
dx      =   L/(nx-1);     % Spacing of grid
x       =   -L/2:dx:L/2; % Grid

% Setup initial temperature profile
T       =   ones(size(x))*Trock;
T(abs(x)<=W/2) = Tmagma;

% set boundary values
T(x == -L/2) = Tleft;
T(x ==  L/2) = Tright;

figure(1); clf;
figure(2); clf;

time = 0;
Trecord = zeros(nt,1);
for n=1:nt               % Timestep loop

    s = kappa * dt / dx^2;   % Heat equation
    
    % Coefficient matrix A
    A = sparse(nx,nx);
    for i = 2:nx-1
        A(i, i-1)  = -s;
        A(i, i  )  = (1+2*s);
        A(i, i+1)  = -s;
    end
    
    % boundary conditions
    A(1,1)   = 1;
    A(nx,nx) = 1;
    
    rhs = zeros(nx,1);
    rhs(2:nx-1) = T(2:nx-1);
    rhs(1)      = Tleft;
    rhs(nx)     = Tright;
    
    % solve linear system of equations
    T = A\rhs;

    % update time
    time  =  time+dt;
    
    % Record T at 5 m distance
    Trecord(n) = T(x==5);
    
    if ~mod(n,10)
        %Plot solution
        figure(1), clf
        plot(x,T);
        xlabel('x [m]')
        ylabel('Temperature [^oC]')
        title(['Temperature evolution after ',num2str(time/day),'days'])
        drawnow
        
        figure(2);
        plot(time/day,Trecord(n),'ro'); hold on;
    end
    
end

