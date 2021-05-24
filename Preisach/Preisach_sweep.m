%% Run simulation
function Out = Preisach_sweep(Grids,In,n)

steps = In.simsteps;
flipped = zeros(steps,1);
flipindex = zeros(Grids.nclusters,3);

if n==0
    Grids.grid = ones(Grids.nnsize,1);
    Out.grids = zeros(steps,Grids.nnsize);
    fields = In.fields;
else
    Grids.grid = Grids.grids(n,:)';
    fields = -In.fields;
end

counter = 0;
for i=1:steps
    In.Eappl = fields(i);
    Grids = induced_energies_cluster(Grids,In);
    dU = Grids.dU;
    flips = find(dU>0);
    if(~isempty(flips))
        flips = setdiff(flips,flipindex(:,1)); %Prevent backflips
        for j=1:length(flips)
            indices = Grids.clusters{flips(j)};
            Grids.grid(indices) = -Grids.grid(indices);
            flipped(i) = flipped(i)+sum(Grids.grid(indices));
            counter = counter+1;
            flipindex(counter,:) = [flips(j) i sum(Grids.grid(indices))];
        end
    end
    if n==0
        Out.grids(i,:) = Grids.grid';
    end  
end

Out.grid = Grids.grid;
Out.flipindex = flipindex;
Out.flipped = flipped;

end
