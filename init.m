function Grids = init(In)
% INIT initializes the Grids variable which contains all dipole configurations
% Grids = INIT(In)

Grids = struct();

%% Sets the initial state and alphas
Grids.nnsize = In.nx*In.ny*In.nz*3;
Grids.grid= ones(Grids.nnsize,1);
Grids.Alpha = permute(repmat(In.alpha*(0:In.nz-1)',1,In.ny,In.nx),[3 2 1]);   %Lefthanded helix with angle In.alpha

%% Draw randoms from parabolic (energy minimum) distribution
x = linspace(-1,1,1000);
probs = -x.^2+1;
[cdf, mask] = unique(cumsum(probs/sum(probs))); %Normalize, integrate and remove doubles
randoms = In.sigmaxy*interp1(cdf,x(mask),rand(In.nx*In.ny,3)); %Invert cdf with random input
totaloffset = repmat(randoms,In.nz*3,1); %random displacement for every column

%     x1 = random('Normal',0,1,In.nx*In.ny,3); % Get randoms from normal distribution

%% Form subcolumns/clusters
alphaoffset = zeros(In.nx,In.ny,In.nz);
randomxy = In.sigmaxy*interp1(cdf,x(mask),rand(In.nx*In.ny,In.nz,2)); %Invert cdf with random input, this matrix is larger then necessary

randomz = In.sigmaz*(1+1/10*rand(In.nx*In.ny,In.nz,1));
randoms = cat(3,randomxy,randomz);

randomlength = In.collength+round(In.sigma*randn(In.nx*In.ny,In.nz)); %Generate subcolumn lengths, this matrix is larger then necessary
randomlength(randomlength<1) = 1; %Set negative column lengths to 1
nmol = In.nx*In.ny*In.nz; %Amount of molecules
Grids.clusters = cell(Grids.nnsize,1); %Contains cluster info, this cell array is larger then necessary
Grids.reversecluster = zeros(Grids.nnsize/3,1); %Contains the number of the cluster a dipole belongs to
Grids.chirality = zeros(In.nx,In.ny,In.nz); %Contains the chirality of each subcolumn
clusterind = 1;

for i=1:In.nx*In.ny %The bottom plane will always be the start of a column
    count = 1;
    clustercount = 1;
    ind = 1; % Random number index
    offset = squeeze(randoms(i,1,:)).'; %start xy-offset
    aoffset = 360*rand(); % Start angle offset
    
    if rand()<0.5 % Set random chirality
        lefthanded = 0;
    else
        lefthanded = 1;
    end
    
    for j=1:In.nz
        if count <= randomlength(i,ind) % Still in subcolumn
            count = count+1;
        else % Column end, go to the next subcolumn with a new offset
            ind = ind + 1;
            count = 1;
            clusterind = clusterind + 1; %Also start new cluster
            clustercount = 1;
            zoffset = offset(3);
            offset = squeeze(randoms(i,j,:)).'; % Set random xy offset
            offset(3) = offset(3)+zoffset;
            aoffset = 2*pi*rand(); % Set random angle offset
            if rand()<0.5 % Set random chirality
                lefthanded = 0;
            else
                lefthanded = 1;
            end
        end
        
        clustercount = clustercount + 1;
        
        Grids.clusters{clusterind} = [Grids.clusters{clusterind} i+(j-1)*In.nx*In.ny];
        Grids.reversecluster(i+(j-1)*In.nx*In.ny) = clusterind;
        
        % Save all offsets in the correct matrices
        [xi,yi] = ind2sub([In.nx In.ny],i);
        alphaoffset(xi,yi,j) = alphaoffset(xi,yi,j) + aoffset;
        Grids.chirality(xi,yi,j) = lefthanded;
        column = i+In.nx*In.ny*(j-1);
        totaloffset(column:nmol:(column+2*nmol),:,1) = totaloffset(column:nmol:(column+2*nmol),:,1) + repmat(offset,3,1); %offset all three dipoles in the molecule
    end
    clusterind = clusterind + 1;
end

Grids.clusters = Grids.clusters(~cellfun('isempty',Grids.clusters));  %Remove empty cells at end of array
Grids.nclusters = length(Grids.clusters);
for i = 1:Grids.nclusters
    Grids.clusters{i} = [Grids.clusters{i}'; Grids.clusters{i}'+nmol ; Grids.clusters{i}'+2*nmol ;]; %Go from molecule index to dipole index
end
Grids.reversecluster = repmat(Grids.reversecluster,3,1); %Go from molecule index to dipole index

Grids.Alpha = Grids.Alpha + alphaoffset; % Apply alpha offset

%% Create Posgrid, the position and dipole grid
Grids.posgrid = zeros(Grids.nnsize,3,2);
for i=1:Grids.nnsize
    [x1,x2,x3,x4] = ind2sub([In.nx In.ny In.nz 3],i);
    Grids.posgrid(i,:,:) = position([x1 x2 x3 x4],Grids,In); % The position function contains the geometric information
end
Grids.posgrid(:,:,1) =  Grids.posgrid(:,:,1) + totaloffset; % Apply offset

%% Replicate sim area
nnsize = Grids.nnsize;
original = Grids.posgrid(:,:,1); % Original box
ind = 1;
periods = zeros(size(original,1)*26,3); % Periodic replications
for i=-1:1
    for j=-1:1
        for k=-1:1
            if ~(i==0&&j==0&&k==0)
                periods(ind:ind+nnsize-1,:) = [original(:,1)-i*In.a*In.nx-j*In.a*In.nx/2 original(:,2)-j*sqrt(3)/2*In.a*In.ny  original(:,3)-k*max(Grids.posgrid(1:Grids.nnsize,3,1))];
                ind = ind+nnsize;
            end
        end
    end
end
periods = [original;periods]; % Add perioidic replications to original box

%% (reverse) Neighbour search
fprintf('Doing (reverse) neighbour search, ')
[IDX,~] = knnsearch(periods,original,'k',In.range+1); % Find nearest neighbours, range+1 because self is also counted
neighbourindexfull = IDX(:,2:end); % Remove self count
Grids.neighbourindex = mod(neighbourindexfull-1,nnsize)+1; % Apply periodic boundaries

Grids.reverseneighbour = cell(nnsize,1);
for i=1:nnsize
    for j=1:In.range
        Grids.reverseneighbour{Grids.neighbourindex(i,j)} = [Grids.reverseneighbour{Grids.neighbourindex(i,j)}; [i,j]];
    end
end

%% Cluster Neighbour search
Grids.neighbourclusters = cell(Grids.nclusters,1);
for i=1:Grids.nclusters
    temp = cell2mat(Grids.reverseneighbour(Grids.clusters{i},:));
    temp = Grids.reversecluster(temp(:,1));
    temp2 = Grids.reversecluster(Grids.neighbourindex(Grids.clusters{i},:));
    temp3 = [temp ; temp2(:)];
    temp3(temp3==i)=[];
    Grids.neighbourclusters{i} = unique(temp3);
end

%% Calculate neighbour distances
fprintf('Calculating interactions, \n')
Grids.neighbourpos = zeros(nnsize,In.range,3); %Saved in Grids for debugging purposes
Grids.dipoles = zeros(nnsize,In.range,3);
for i=1:nnsize
    Grids.neighbourpos(i,:,:) = periods(neighbourindexfull(i,:),:); % Positions of neighbours
    Grids.dipoles(i,:,:) = Grids.posgrid(Grids.neighbourindex(i,:),:,2); % Dipoles of neighbours
end

pos = permute(repmat(Grids.posgrid(:,:,1),1,1,In.range),[1 3 2]); % Get all the position vectors in the right format

D = Grids.neighbourpos - pos; % Distance vector
Grids.r = 1./sqrt(sum(D.^2,3)); % Inverse distance [nm^-3]
Grids.ru = bsxfun(@times, D, Grids.r); % Unit distance vector
Grids.r = (Grids.r*1e9).^3; % Inverse distance in correct units [m^-3]

end