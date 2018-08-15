%{
Model that simulates the flipping of dipoles during a Preisach measurement
of a BTA ferroelectric.
Tim Cornelissen, Linkoping University, Summer 2018

                    YAao,
                    Y8888b,
                  ,oA8888888b,
            ,aaad8888888888888888bo,
         ,d888888888888888888888888888b,
       ,888888888888888888888888888888888b,
      d8888888888888888888888888888888888888,
     d888888888888888888888888888888888888888b
    d888888P'                    `Y888888888888,
    88888P'                    Ybaaaa8888888888l
   a8888'                      `Y8888P' `V888888
 d8888888a                                `Y8888
AY/'' `\Y8b                                 ``Y8b
Y'      `YP                                    ~~
         `'
%}

close all
rng('shuffle');
finishup = onCleanup(@() delete(findall(0,'type','figure','tag','TMWWaitbar')));

%% Input Parameters
% Lattice parameters
In.nx = 8; %[#] molecules in x direction
In.ny = 8; %[#] molecules in y direction
In.nz = 30; %[#] molecules in z direction
In.alpha = pi/3; %[radian] The rotation of the dipoles configuration between z-planes, pi/3
In.a = 3; %[nm] Hexagonal packing distance, C6=1.67, C10=2, C18=2.55
In.L = 0.28; %[nm] Distance from dipole to center of molecule ~0.28nm
In.c = 0.35;  %[nm] Intercolumnar packing distance, ~0.35nm
In.beta = 2*pi/9; %[radian] Out of plane rotation of dipoles
In.range = 40; %[#] Interaction range
In.collength = 7; %[#] average subcolumn/cluster length
In.sigma = 3; %[#] disorder in column length
In.sigmaxy = 0.02; %[nm] disorder in xy direction
In.sigmaz = 0.02; %[nm] disorder in z direction

% Material and environment parameters
In.mu = 4; %[D] Dipole strength, ~4D
In.eps = 2; %[-] Relative dielectric constant 
In.Eappl = 0; %[V/m] Applied electric field
In.T = 450; %[K] Temperature
In.pol = 1; %[eA^2/V] Polarizability

% Simulation parameters
In.filename = '/filename.mat'; % Filename
In.simsteps = 10; % Steps for simulation
In.iter = 3; % Amount of iterations to solve induced dipole interactions, 3-4 is sufficient
In.fields = linspace(0,-2e9,In.simsteps); %[V/m] Field sweep

%% Conversion of units and physical constants
In.eps0 = 8.854187816e-12; %[F/m] Vacuum permittivity
In.ke = 8.987551787368176e9/In.eps; %[N m^2 C^-2] Coulombs constant 1/(4pi*eps0)
In.pol = In.pol/6.24e18/1e20; %[Cm^2/V] Polarizability
In.mu = In.mu*3.33564e-30; %[C m]
In.B = 1/(In.T*1.380650324e-23); %[1/J] 1/kT
In.fieldconstant = -In.B*In.mu^2*In.ke; % [m^3] For faster field calculations, Divide by eps or not??!?
In.reaction = (In.eps-1)/(2*In.eps-1); % [-] Reaction field constant
In.filename = [pwd In.filename];

%% Check for version
if verLessThan('matlab', '9.1')
    return % Code doesn't work on <2016b due to use of implicit expansion i.s.o. bsxfun
end

%% First sweep
fprintf('Starting first sweep')
Grids = init(In); %initiliaze grid
Out = Preisach_sweep(Grids,In,0);

%% Run all sweeps back in parallel
Grids.grids = Out.grids;
fprintf('Starting sweeps back')
fprintf(1,'%s\n\n',repmat('.',1,In.simsteps)); %Progress indication

parfor i=1:In.simsteps 
    Output(i) = Preisach_sweep(Grids,In,i);
    fprintf(1,'\b.\n'); 
end

%% Process and show results
% Show morphology
figure
scatter3(Grids.posgrid(1:Grids.nnsize,1,1),Grids.posgrid(1:Grids.nnsize,2,1),Grids.posgrid(1:Grids.nnsize,3,1),[],[-Grids.grid -Grids.grid zeros(Grids.nnsize,1)]);
axis equal

% Cumulative Preisach distribution
data = zeros(In.simsteps );
for i = 1:In.simsteps 
    data(i,:) = Output(i).flipped';
end
cumP = cumsum(data,2);

[X,Y] = meshgrid(-In.fields,In.fields);

figure
surf(X,Y,cumP)

%% Saving
save(In.filename,'In','Grids','Output')
fprintf('Result saved')
