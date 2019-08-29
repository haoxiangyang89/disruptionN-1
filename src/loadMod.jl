# load functions and packages

using JuMP, Ipopt, Gurobi, Combinatorics, JLD, HDF5, DelimitedFiles;

include("def.jl");
include("readin.jl");
include("main.jl");
include("forwardPass.jl");
include("backwardPass.jl");
