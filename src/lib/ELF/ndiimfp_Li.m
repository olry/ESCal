function [diimfp,dsep,reduced_diimfp] = ndiimfp_Li(osc,E0,depth,alpha,decdigs,varargin)

%%
%{
   Calculates the normalised DSEP (in eV^-1*A^-1) and DIIMFP (in eV^-1)
   for a given energy, angle and depth
   from solid to vacuum
   according to the Li algorithm Eq.(9) from
   Y.C. Li et al. / Surface Science 589 (2005) 67-76.
%}

if nargin<5, decdigs=10; end
if nargin<4
    warning ('Error in input format')
else   
    %% Clear bulk
    if depth >= 0
        x_in_b = ndiimfp(osc,E0,10,false);
    else
        x_in_b = zeros(size(osc.eloss));
    end

    %% Reduced bulk

    theta = 0:pi/2/5:pi/2;
    phi = 0:2*pi/10:2*pi;
    v = sqrt(2*E0/h2ev);

    q = zeros(length(osc.eloss),2^(decdigs-1)+1);
    Im = zeros(length(osc.eloss),2^(decdigs-1)+1,length(theta));
    x_in_s = zeros(size(osc.eloss));

    q_minus = sqrt(2*E0/h2ev) - sqrt(2*(E0/h2ev-osc.eloss/h2ev));
    q_plus =  sqrt(2*E0/h2ev) + sqrt(2*(E0/h2ev-osc.eloss/h2ev));

    for i = 1:2^(decdigs-1)+1
        q(:,i) = q_minus + (i-1)*(q_plus-q_minus)/2.^(decdigs-1);
    end   
    if strcmp(osc.model,'Mermin')
        q(q==0) = 0.01;
    end

    Q = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(sin(theta),1,1,[]));
    indQ = bsxfun(@eq,Q,0);
    Q(indQ) = eps;

    qz = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(cos(theta),1,1,[]));

    v_perpendicular = cosd(alpha).* v;

    r = depth./a0; %./cosd(alpha);

    qz_r_cosalpha = qz.*cosd(alpha).*r;
    q_sinsquared_theta = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(sin(theta).^2,1,1,[]));

    top = q_sinsquared_theta .* cos(qz_r_cosalpha) .* exp(-abs(r).*Q.*cosd(alpha));

    % w_wave
    qv_sintheta = bsxfun(@times,repmat(q .* v, 1, 1, length(theta)),reshape(sin(theta),1,1,[]));
    exdim = repmat(qv_sintheta, 1, 1, 1, length(phi)); %add extra dimension over phi
    qv_sintheta_cos_phi = bsxfun(@times,exdim,reshape(cos(phi),1,1,1,[]));    
    w_wave = bsxfun(@minus,repmat((osc.eloss'/h2ev)',1,2^(decdigs-1)+1,length(theta),length(phi)), qv_sintheta_cos_phi.*sind(alpha));

    Qv_perpendicular_squared = Q.^2.*v_perpendicular^2;

    bottom = bsxfun(@plus,w_wave.^2,Qv_perpendicular_squared);

    %      exdimdep = repmat(Q, 1, 1, 1, length(phi),length(r)); %add extra dimension over r
    %      exprQ = bsxfun(@times,reshape((-1)*abs(r),1,1,1,1,[]),exdimdep).*cosd(alpha);
    %      top_in = bsxfun(@times,bsxfun(@times,q_sinsquared_theta,cos(qz_r_cosalpha)),exp(exprQ));

    for i = 1:length(theta)
        osc.qtran = Q(:,:,i)/a0;
        Im(:,:,i) = eps_sum_allwq(osc,'bulk');
    end

    %      topbot = bsxfun(@rdivide,top_in,repmat(bottom, 1, 1, 1, 1, length(r)));
    topbot = bsxfun(@rdivide,top,bottom);
    if depth >= 0 
        romall_in = bsxfun(@times,Im,topbot);
        romall_in(isnan(romall_in)) = 0;
        romall_in = romall_in.*stepfunction(r);
        res_in = squeeze(trapz(theta,trapz(phi,romall_in,4),3));

        x_in_b_reduced = zeros(size(osc.eloss));

        for i = 1:length(osc.eloss)
            x_in_b_reduced(i) = (-2)*cosd(alpha)/(pi^3) * trapz(q(i,:),res_in(i,:),2)*1/h2ev/a0;
        end
    else
        x_in_b_reduced = zeros(size(osc.eloss));
    end
     
    %% Surface
    
    for i = 1:length(theta)
        osc.qtran = Q(:,:,i)/a0;
        Im(:,:,i) = eps_sum_allwq(osc,'surface');
    end
   
    %================= inside ===================
    romall_in = bsxfun(@times,Im,topbot);
    romall_in(isnan(romall_in)) = 0;
    romall_in = romall_in.*stepfunction(r);
    res_in = squeeze(trapz(theta,trapz(phi,romall_in,4),3));

    %================= outside ==================
    top = q_sinsquared_theta .* exp(-abs(r).*Q.*cosd(alpha));
    topbot = bsxfun(@rdivide,top,bottom);
    
    add = 2.*cos(w_wave.*r./v) - exp(-abs(r).*Q.*cosd(alpha));

    romall_out = bsxfun(@times,Im,topbot).*add;
    romall_out(isnan(romall_out)) = 0;
    romall_out = romall_out.*stepfunction(-r);
    res_out = squeeze(trapz(theta,trapz(phi,romall_out,4),3));

    for i = 1:length(osc.eloss)
        x_in_s(i) = 4*cosd(alpha)/(pi^3) * (trapz(q(i,:),res_in(i,:)) + trapz(q(i,:),res_out(i,:)))*1/h2ev/a0;
    end   

    
    %%     
    dsep = x_in_s; % only surface component
    diimfp = x_in_b; % clear bulk  
    reduced_diimfp = x_in_b_reduced; % reduced bulk
end

    %% Heaviside function
    function x = stepfunction(depth)
        if depth > 0
            x = 1;
        elseif depth < 0
            x = 0;
        elseif depth == 0
            x = 0.5;
        end
    end

end










