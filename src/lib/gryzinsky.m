function sigma = gryzinsky(E,Eb,shell)
    n_ele = [2 2 2 4 2 2 4 4 6 2 2 4 4 6 6 8 2 2 4 4 6 2 2 4];
    u = E/Eb;
    e2 = 14.39988; % square of elem. charge in eV Angstroem
    sigma = n_ele(shell)*e2*e2*pi*1/u*((u-1)/(u+1)^1.5)...
        *log(2.7+sqrt(u-1))*(1+2/3*(1-0.5/u))/Eb/Eb;
end