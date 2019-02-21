function sigma = bote(E,Eb,a1,a2,a3,a4,a5)
    u = E/Eb;
    sigma = 4*pi*a0^2*((u-1)/u^2)*(a1 + a2*u + a3/(u+1) + a4/((u+1)^3) + a5/((u+1)^5))^2;
end