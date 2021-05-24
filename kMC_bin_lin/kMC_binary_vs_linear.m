% A quick test for a kMC simulation of a simple Ising-like model
% In kMC the limiting factor is searching for the cell to flip
% We compare a linear search to a binary tree search approach
% The kMC algorithm is nonsense, but a MWE to make the efficiency comparison
% Tim Cornelissen, July 2020

%% Initialize
close all % Close all figures

Grid_size = 5000; % Size of simulation grid
Rates = round(rand(1,Grid_size)); % Flipping rates are just a random boolean
Range = 30; % Interaction range (# neighbors)
N_steps = 1e4; % # Simulation Steps
N_sims = 100; % Number of trial simulations
Error = false; % Error switch

%% Determine neighbors
Neighbors = zeros(Grid_size,Range); % Initialize matrix which will contain who are a cell's neighbors
for k = 1:Grid_size
    Neighbors(k,:) = (k+1):(k+Range); % Each cell has Range neighbors on its right
end
Neighbors = mod(Neighbors-1,Grid_size)+1; % Apply periodic boundaries

%% Build Binary tree
Tree_height =  ceil(log2(Grid_size)); % Height of the binary tree
Rates_tree = zeros(Tree_height+1,Grid_size+1); % Initialize binary tree with flipping rates
Rates_tree(end,1:end-1) = Rates; % Fill lowest level of tree with flipping rates

for k = Tree_height:-1:1 % Step through levels of tree
    for j = 1:ceil(Grid_size/(2*(Tree_height+1-k))) % Tree halves in width each level
        Rates_tree(k,j) = Rates_tree(k+1,2*j-1) + Rates_tree(k+1,2*j); % Node value is sum of children nodes
    end
end
Rates_tree(:,end) = []; % Clear row of zeros at the end

%% Run and compare sims
Rands = rand(N_steps,N_sims); % Generate sets of random numbers for simulation
Simtime_lin = zeros(N_sims,1); % Initialize simtime array
Simtime_bin = zeros(N_sims,1); % Initialize simtime array

for k = 1:N_sims
    tic
    flipped_lin = kMC_lin(N_steps,Rands(:,k),Rates,Neighbors); % Linear search sim
    Simtime_lin(k) = toc; % Save simulation time
    
    tic
    flipped_bin = kMC_bin(N_steps,Rands(:,k),Rates_tree,Neighbors,Tree_height); % Binary search sim
    Simtime_bin(k) = toc; % Save simulation time
    
    if sum(abs(flipped_lin-flipped_bin)) > 0 % Check if results are the same
        Error = true;
        break;
    end
end

%% Show results
if(~Error) % Make histogram of simulation times
    figure
    histogram(Simtime_lin)
    hold on
    histogram(Simtime_bin)
    legend('Linear', 'Binary')
    xlabel('Simulation time (s)')
    ylabel('Counts')
else
    fprintf('Error: algorithms give different results')
end

