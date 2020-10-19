function dsep = dsep_Tung_2(osc,E0,decdigs,norm,theta,varargin)
%%
%{
   Calculates the normalised DSEP
   for a given energy and angle
   according to the Tung algorithm Eq.(9)
   C.J. Tung et al. / Phys. Rev. B 49 (1994) 16684.
%}

if nargin<4, decdigs=10; end
if nargin<3
    warning ('Error in input format')
else
    C = 137.036;
    cur_rel_cor_factor = (1 + E0/h2ev/ (C*C))*(1 + E0/h2ev / (C*C)) / (1 + E0/h2ev/ (2.0*C*C));
    ind = bsxfun(@gt,osc.eloss,0);
    eloss = osc.eloss; 
    
    q1 = sqrt(E0/h2ev*(2 + E0/h2ev / (C*C))) - sqrt((E0/h2ev - osc.eloss(ind)/h2ev).*(2.0 + (E0/h2ev - osc.eloss(ind)/h2ev) / (C*C)));  % relativistic integration limit, see Shinotsuka SIA 47 871
    q2 = sqrt(E0/h2ev*(2 + E0/h2ev / (C*C))) + sqrt((E0/h2ev - osc.eloss(ind)/h2ev).*(2.0 + (E0/h2ev - osc.eloss(ind)/h2ev) / (C*C)));
    
    q = zeros(length(osc.eloss(ind)),2^(decdigs-1)+1);
    
    for i = 1:2^(decdigs-1)+1
        q(:,i) = q1 + (i-1)*(q2-q1)/2.^(decdigs-1);
    end
    
    if strcmp( osc.model,'Mermin')
        q(q==0) = 0.01;
    end
    
    osc.qtran = q/a0;
    osc.eloss = eloss(ind);
    osc.qtran = 0.01;
    epsilon = epsilon_allwq(osc);
    epsfraction = (epsilon - complex(1, 0)).*(epsilon - complex(1, 0)) ./ (epsilon.*(epsilon + complex(1, 0)));
    qs = sqrt(abs(q.^2 - (osc.eloss/h2ev + 0.5*q.^2).*(osc.eloss/h2ev + 0.5*q.^2) / (2.0 * E0/h2ev)));    
    result = imag(epsfraction) .*qs ./ (q.^3);
    
    rel_cor_factor = (1 + E0/h2ev / (C*C))*(1 + E0/h2ev / (C*C)) / (1 + E0/h2ev / (2.0*C*C));
    
    tointegrate = 2*rel_cor_factor*result / (pi*E0/h2ev);
    
    dsep = zeros(size(osc.eloss));
    
    for i=1:length(osc.eloss)
        dsep(i) = cur_rel_cor_factor * trapz(q(i,:),tointegrate(i,:)) /h2ev;
    end
    dsep(1) = eps;
    ind = osc.eloss < 1;
    dsep(ind) = nan;
    dsep(1) = eps;
    fill = fillmissing(dsep,'spline');
    dsep = fill;
    if norm
        dsep = dsep / trapz(osc.eloss,dsep);
    end
end

end

