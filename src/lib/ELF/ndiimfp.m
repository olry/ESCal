function diimfp = ndiimfp(osc,E0,decdigs,norm,gapEffect,relativistic,varargin)
%%
%{
   Calculates the normalised DIIMFP
   for a given energy.
%}
if nargin<6, relativistic = false; end
if nargin<5, gapEffect = false; end
if nargin<4, norm = true; end

if nargin<3, decdigs = 10; end
if nargin<2
    warning ('Error in input format')
else
    E0 = E0/h2ev;
    omega = osc.eloss / h2ev;
    if relativistic
        C = 137.036; % a.u.
        q1 = sqrt(E0*(2 + E0 / C^2 )) - sqrt((E0 - omega).*(2.0 + (E0 - omega) / C^2 ));  % relativistic integration limit, see Shinotsuka SIA 47 871
        q2 = sqrt(E0*(2 + E0 / C^2 )) + sqrt((E0 - omega).*(2.0 + (E0 - omega) / C^2 ));
        q = zeros(length(omega),2^(decdigs-1)+1);
    
        for i = 1:2^(decdigs-1)+1
            q(:,i) = q1 + (i-1)*(q2-q1)/2.^(decdigs-1);
        end

        if strcmp( osc.model,'Mermin')
            q(q==0) = 0.01;
        end

        x_in = zeros(size(osc.eloss));

        osc.qtran = q/a0; 
        epsilon = epsilon_allwq(osc,gapEffect);
        res = imag(-1./epsilon)./q;
        rel_cor_factor = (1 + E0/C^2)^2 / (1 + E0/(2*C^2));

        x_in(1) = eps;
        for i = 2:length(osc.eloss)
            x_in(i) = trapz(q(i,:),res(i,:));
        end

        x_in = rel_cor_factor*x_in / (pi*E0) / (h2ev*a0);
    else           
        qmin = log( sqrt(2*E0) - sqrt(2*(E0 - omega)) );
        qmax = log( sqrt(2*E0) + sqrt(2*(E0 - omega)) );
        q = zeros(length(osc.eloss),2^(decdigs-1)+1);
        x_in = zeros(size(osc.eloss));

        for i = 1:2^(decdigs-1)+1
            q(:,i) = qmin + (i-1)*(qmax-qmin)/2.^(decdigs-1);
        end

        if strcmp( osc.model,'Mermin')
            q(q==0) = 0.01;
        end

        osc.qtran = exp(q)/a0; 
        res = eps_sum_allwq(osc,'bulk',gapEffect);

        x_in(1) = eps;
        for i = 2:length(osc.eloss)
            x_in(i) = 1/pi/E0 * trapz(q(i,:),res(i,:))*(1/h2ev/a0);
        end
    end
    if (norm)
        diimfp = x_in./trapz(osc.eloss,x_in); %normalized
    else
        diimfp = x_in; %for imfp calculation
    end
end










