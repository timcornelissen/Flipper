function Grids = induced_energies_cluster(Grids,In)
% INDUCED_ENERGIES calculates the energy difference for flipping the clusters
% Grids = INDUCED_ENERGIES(Grids,In)

% intialize 
U_start = zeros(Grids.nclusters,2);


for i = 1:Grids.nclusters
   
    clusterind = Grids.neighbourclusters{i};
    
    indices = [cell2mat(Grids.clusters(clusterind)); Grids.clusters{i}] ;  % Indices of all involved dipoles
    
    basegrid = zeros(Grids.nnsize,1);
    basegrid(indices) = Grids.grid(indices);
    
    % Make filter to select only neighbours
    ind_temp = Grids.neighbourindex(indices,:);
    ind = zeros(size(ind_temp));
    for j = 1:length(indices)
        ind(ind_temp==indices(j) )= j;
    end
    filter = (ind~=0);
    ind(~filter) = 1;
    ind = cat(3,ind,ind+length(indices),ind+2*length(indices));
    
    %U, cluster i in original state
    grid = basegrid;
    U_start(i,1) = solve_induced(Grids,In,indices,grid,ind,filter);
    
    %U, cluster i flipped
    grid(Grids.clusters{i}) = -grid(Grids.clusters{i});
    U_start(i,2) = solve_induced(Grids,In,indices,grid,ind,filter);
    
end

Grids.dU = -diff(U_start,1,2);

end