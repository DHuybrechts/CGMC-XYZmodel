function [ time, Sx1, Ct, rho, Mx, My, Mz, MX, MY, MZ] = CalculateSssxx( gamma, Jx, Jy, Jz, Nxc, Nyc, Nx, Ny, T, dt)
%-------------------------------------------------------------------------%
%Calculate the time evolution of the quantum trajectory and calculate
%the steady-state spin structure factor.
%Parameters:
%   gamma       system parameter
%   Jx,Jy,Jz    system parameter
%   Nxc         number of sites in x direction of cluster
%   Nyc         number of sites in y direction of cluster
%   Nx          number of rows
%   Ny          number of columns
%   T           Final time of quantum trajectory.
%   dt          Time steps to save.
%-------------------------------------------------------------------------%
    %Some parameters:
    N = Nx*Ny;                                                                                          %number of sites in lattice
    Nc = Nxc*Nyc;                                                                                       %number of sites in cluster
    NC = N/Nc;                                                                                          %number of clusters in lattice
    N_coeff = 2^Nc;                                                                                     %number of coefficients in cluster wave function
    
    Cin = ones(N_coeff*NC,1)/sqrt(N_coeff);                                                             %initial magnetization is 1 in the x direction
  
    if Nxc == 1 && Nx ==1
        [~, ~, out] = CalcClusterNeighbours1D( Nxc, Nyc, Nx, Ny );
    else
        [~, ~, out] = CalcClusterNeighbours( Nxc, Nyc, Nx, Ny );                                        %calculates 'all' neirest neighbour as well as 'in' and 'out' of the cluster.
    end
    clustconfig = ClusterConfiguration(Nxc, Nyc, Nx, Ny);                                               %matrix where the row index represents the cluster and contains the indices of the sites in that respective cluster.
    sig = GetAllOperators(Nc);                                                                          %Calculates all cluster operator acting on the different sites in the cluster.
    sig = nDMatrixToCellArraySparse(sig);                                                               %Transform nD-matrices to cell arrays
    [switch_xy, magn_y, magn_z] = CalcAlternativeOperator(sig);                                         %Uses the permutation representation of the operators and a vector for pointwise multiplication with sign changes etc.
    NNM = NeirestNeighbourMatrix(out);                                                                  %Connection matrix, each row is a site on the lattice and the columns with a 1 are the neirest neighbours outside the cluster.
    A1 = H1( gamma, Jx, Jy, Jz, Nxc, Nyc, Nx, Ny );                                                     %Time independt part of the Hamiltonian
    Tnorm = kron(eye(NC),ones(1,N_coeff));                                                              %Matrix to help calculate the norm of a product wave function.
    
    %data variables
    sx = zeros(floor(T/dt)+1,1);                                                                        %steady-state spin structure factor
    mx = sx;                                                                                            %mean x magnetization of lattice at certain time
    my = sx;                                                                                            %mean y magnetization of lattice at certain time
    mz = sx;                                                                                            %mean z magnetization of lattice at certain time
    MX = zeros(floor(T/dt)+1,N);                                                                        %x magnetization of each lattice site at certain time                    
    MY = MX;                                                                                            %y magnetization of each lattice site at certain time
    MZ = MX;                                                                                            %z magnetization of each lattice site at certain time
    RHO = cell(1, floor(T/dt)+1);                                                                       %Density matrix at each time t
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t_left = 10;                                                                                        %Time evolution untill this time to remove the dependence on the initial wave function.
    check = true;                                                                                       %Check variable to determine if timestep dt endend (Check = false).
    e1 = rand;                                                                                          %Random variable to determine jump time.
    while check
        
        [t_s, Cin, ~, ~, ~] =EffH(gamma, Jx, Jy, Jz, Nxc, Nyc, Nx, Ny,...
                clustconfig, NNM, sig, switch_xy, magn_y, magn_z, A1, Cin, Tnorm, [0 t_left], e1);      %Time evolution up untill t_left or untill the jump time t_j.
        if t_s(end) == t_left                                                                           %NO JUMP in the chosen time interval.
            Cin = transpose(Cin(end,:))./NormC(transpose(Cin(end,:)), NC);                              %Normalize wave function
            check = false;                                                                              %Check becomes false to go to the next time step.
        else                                                                                            %JUMP occurs on time t_j = t_s(end).
            i_jump = GetSite(ChanceInterval(transpose(Cin(end,:)), clustconfig, sig, Tnorm), rand);     %Determine the site that jumps using their probability distribution.
            Cin = Jump(transpose(Cin(end,:)), clustconfig, sig, i_jump)./...
                  NormC(Jump( transpose(Cin(end,:)), clustconfig, sig, i_jump ), NC);                   %Update wave function with jump and normalize.
            t_left = t_left - t_s(end);                                                                 %Time left untill next dt step.
            e1 = rand;                                                                                  %Pick a new random number to determine jump time.
        end
    end
    
    %calculate some expectation values and save them:
    C = transpose(Cin);
	sx(1) = CalcSss(Cin, sig, clustconfig, N);                                                          %Calculate the spin structure factor.
    Sigxyz = CalcExpSig(Cin, switch_xy, magn_y, magn_z, clustconfig);                                   %Calculate the expectation values
    mx(1) = mean(Sigxyz(:,1));                                                                          %mean x magnetization at this time
    my(1) = mean(Sigxyz(:,2));                                                                          %mean y magnetization at this time
    mz(1) = mean(Sigxyz(:,3));                                                                          %mean z magnetization at this time
    MX(1,:) = transpose(Sigxyz(:,1));                                                                   %x magnetization at each site at this time
    MY(1,:) = transpose(Sigxyz(:,2));                                                                   %y magnetization at each site at this time
    MZ(1,:) = transpose(Sigxyz(:,3));                                                                   %z magnetization at each site at this time
    RHO{1} = CalcRho(Cin, Nxc, Nyc, Nx,Ny);                                                             %Calculate the density matrix from the wave function
    
    %Continue time evolution of the trajectory
    for i = 1:floor(T/dt)
        t_left = dt;                                                                                    %Time left untill next dt
        check = true;                                                                                   %Check variable to determine if timestep dt endend (Check = false).
        e1 = rand;                                                                                      %Random variable to determine jump time.
        while check
            [t_s, Cin, ~, ~, ~] = EffH(gamma, Jx, Jy, Jz, Nxc, Nyc, Nx, Ny,...
                clustconfig, NNM, sig, switch_xy, magn_y, magn_z, A1, Cin, Tnorm, [0 t_left], e1);      %Time evolution up untill t_left or untill the jump time t_j.
            if t_s(end) == t_left                                                                       %NO JUMP in the chosen time interval.
                Cin = transpose(Cin(end,:))./NormC(transpose(Cin(end,:)), NC);                          %Normalize wave function
                C(size(C,1)+1,:) = transpose(Cin);                                                      %Save dt value.
                check = false;                                                                          %Check becomes false to go to the next time step.
            else                                                                                        %JUMP occurs on time t_j = t_s(end).
                i_jump = GetSite(ChanceInterval(transpose(Cin(end,:)), clustconfig, sig, Tnorm), rand); %Determine the site that jumps using their probability distribution.
                Cin = Jump(transpose(Cin(end,:)), clustconfig, sig, i_jump)./...
                    NormC(Jump( transpose(Cin(end,:)), clustconfig, sig, i_jump ), NC);                 %Update wave function with jump and normalize.
                t_left = t_left - t_s(end);                                                             %Time left untill next dt step.
                e1 = rand;                                                                              %Pick a new random number to determine jump time.
            end
        end
        Sigxyz = CalcExpSig(Cin, switch_xy, magn_y, magn_z, clustconfig);
        
        sx(i+1) = CalcSss(transpose(C(i+1,:)), sig, clustconfig, N);                                   %Calculate the spin structure factor.
%         mx(i+1) = CalcSss2(transpose(C(i+1,:)), sig, clustconfig, N);
        mx(i+1) = mean(Sigxyz(:,1));
        my(i+1) = mean(Sigxyz(:,2));
        mz(i+1) = mean(Sigxyz(:,3));
        MX(i+1,:) = transpose(Sigxyz(:,1));
        MY(i+1,:) = transpose(Sigxyz(:,2));
        MZ(i+1,:) = transpose(Sigxyz(:,3));
        RHO{i+1} = CalcRho(transpose(C(i+1,:)), Nxc, Nyc, Nx, Ny);
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    time = (0:dt:floor(T/dt)*dt)';              %time interval
    Sx1 = sx;                                   %steady-state spin structure factor
    Mx = mx;                                    %mean x magnetization
    My = my;                                    %mean y mangetization
    Mz = mz;                                    %mean z magnetization
    Ct = C;                                     %wave function
    rho = RHO;                                  %density matrix
end