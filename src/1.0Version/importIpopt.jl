function importIpopt()
    mp = Model(solver = IpoptSolver(linear_solver = "ma27"));
    @variable(mp, x1 >= 0);
    @variable(mp, x2 >= 0);
    cons1 = @constraint(mp, x1 + 2x2 <= 1);
    cons2 = @constraint(mp, x2 + 2x1 <= 1);
    @objective(mp, Max, x1 + x2);
    solve(mp);
end
