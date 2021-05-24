This program does a simple check to compare two different algorithms for the most time-consuming step in a kinetic Monte Carlo (kMC) simulation: finding the event to be carried out. This was helpful for me to determine the most efficient way to program a big kMC simulation for my PhD research.

## The Ising Model
This simple simulation is based on an Ising model. This model considers a 1D chain of cells that contain a dipole (also called spin) that points either up or down. A cell interacts with its neighbors, which gives an energy difference between its up and down state. This energy difference determines the probability that a cell will flip its dipole. Each cell thus has a certain flipping probability, or rate, depending on the configuration of its neighbors. 

## kMC algorithm
To calculate how the configuration of dipoles evolves over time, a kMC simulation can be used:
    Choose an initial state for all N cells
    Calculate the flipping rates r<sub>i of all cells 
    Calculate the cumulative flipping rate array R<sub>j</sub>= &Sigma;<sub>i</sub><sup>j</sup> r<sub>i</sub> for j=1,...,N
    For step = 1:N_steps
        Get a uniform random number u
        Find the cell to be flipped by finding for which cell k we have R<sub>k-1</sub><uR<sub>N</sub><R<sub>k</sub>
        Flip cell k
        Update r<sub>k</sub> and consequently R<sub>j</sub>
    Finish up simulation and calculate properties of final state (such as correlation functions)

The most computationally expensive step is finding the cell to flip.

## kMC_binary_vs_linear
In kMC_binary_vs_linear we compare two approaches to do this search: a ‘linear’ search using Matlab’s find function versus a binary tree search. The binary tree approach require more preparation and bookkeeping, but is faster in the end. 
In this simple example, the underlying Ising-like model is physical nonsense. A cell has  neighbors to its right only, and the flipping rates are simply Booleans. If a cell flips, its flipping rate as well as that of all its neighbors reverses. While this doesn’t make sense for any real physical model, it is sufficiently similar to a real Ising kMC simulation to make the efficiency comparison.
