function Flipped = kMC_bin(N_steps,Rands,Rates_tree,Neighbors,Tree_height)
% Performs a simple kMC simulation with binary tree search, returns array
% of flipped indices
% Flipped = kMC_bin(N_steps,Rands,Rates_tree,Neighbors,Tree_height)

Flipped = zeros(N_steps,1); % Initialize array of flipped indices

for i=1:N_steps
    r = Rands(i)*Rates_tree(1,1); % Random number times total rate
    
    % Walk down tree to find cell to flip
    index = 1; % Start at the left of the tree
    level = 1; % Start at the top level
    
    while level<Tree_height+1
        rate = Rates_tree(level+1,2*index-1); % Rate of node on left branch
        if r <= rate
            index = 2*index-1; % Step down left branch
        else
            r = r - rate; % Subtract rate of node on left branch
            index = 2*index; % Step down right branch
        end
        level = level+1; % Go down one level       
    end
    
    Flipped(i) = index; % Save index of flipped cell
    changed_indices = [index Neighbors(index,:)]; % Rates of flipped cell and neighbors need to be updated
    
    % Walk up tree to update rates
    for l = 1:length(changed_indices)
        level = Tree_height+1; % Start at bottom level
        changed_index = changed_indices(l); % Index of cell rate to be updated
        ratechange = ~Rates_tree(end,changed_index)-Rates_tree(end,changed_index); % Change in rate
        
        while level>0 % Walk up tree to propogate updated rate
            Rates_tree(level, changed_index) = Rates_tree(level,changed_index) + ratechange; % Update rate
            level = level-1; % Go up one level
            changed_index = ceil(changed_index/2); % Step to the correct index in the next level
        end
    end
end

end