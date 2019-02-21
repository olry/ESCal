function sigma = casnati(E,Eb,shell)
    n_ele = [2 2 2 4 2 2 4 4 6 2 2 4 4 6 6 8 2 2 4 4 6 2 2 4];
    d0 = -0.0318;
    d1 = 0.3160;
    d2 = -0.1135;
    b0 = 10.57;
    b1 = -1.736;
    b2 = 0.317;
    u = E/Eb;
    Rydberg = 13.605693009; %eV
    exx = d0+d1/u+d2/u/u;
    fi = b0*exp(b1/u+b2/u/u);
    psi = (Eb/Rydberg)^exx;

    sigma = n_ele(shell)*a0*a0*((Rydberg/Eb)^2)*psi*fi*log(u)/u;
end