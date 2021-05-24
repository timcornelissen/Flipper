function Flipped = kMC_lin(N_steps,Rands,Rates,Neighbors)
% Performs a simple kMC simulation with linear search, returns array of 
% flipped indices 
% Flipped = kMC_lin(N_steps,Rands,Rates,Neighbors)

Flipped = zeros(N_steps,1); % Initialize array of flipped indices 

for i = 1:N_steps
    Rate_sum = cumsum(Rates); % Cumulative flipping rates
    index = find((Rands(i)*Rate_sum(end)) < Rate_sum,1); % Linear search for index to flip
    
    % Flip event
    Flipped(i) = index; % Save index of flipped cell
    Rates(index) = ~Rates(index); % Update flip rate of flipped cell
    Rates(Neighbors(index,:)) = ~Rates(Neighbors(index,:)); % Update flip rate of neighbors
end

end