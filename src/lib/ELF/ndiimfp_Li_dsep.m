function [dsep,siimfp] = ndiimfp_Li_dsep(osc,E0,depth,alpha,decdigs,varargin)

%%
%{
   Calculates the normalised DSEP (in eV^-1*A^-1)
   for a given energy, angle and depth
   from solid to vacuum
   according to the Li algorithm Eq.(9) from
   Y.C. Li et al. / Surface Science 589 (2005) 67-76.
%}

if nargin<5, decdigs=10; end
if nargin<4
    warning ('Error in input format')
else
    
    qmin = sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-osc.eloss/h2ev));
    qmax = sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-osc.eloss/h2ev));
    
    q = zeros(length(osc.eloss),2^(decdigs-1)+1);

    for i = 1:2^(decdigs-1)+1
        q(:,i) = qmin + (i-1)*(qmax-qmin)/2.^(decdigs-1);
    end
    
    theta = 0:pi/2/5:pi/2;
    phi = 0:2*pi/10:2*pi;
    
    Im = zeros(length(osc.eloss),2^(decdigs-1)+1,length(theta));
    
    qz = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(cos(theta),1,1,[]));
    
    Q = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(sin(theta),1,1,[]));
    indQ = bsxfun(@eq,Q,0);
    Q(indQ) = 0.0001;
    
    v_per = cosd(alpha).*sqrt(2*E0/h2ev);
    r = depth./a0./cosd(alpha);
    exdimr = repmat(qz, 1, 1, 1,length(phi), length(r)); %add extra dimension over r and phi
    qzrcos = bsxfun(@times,exdimr,reshape(r,1,1,1,1,[])).*cosd(alpha);
    
    qsintheta = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(sin(theta).^2,1,1,[]));
    
    qvsintheta = bsxfun(@times,repmat(q.*sqrt(2*E0/h2ev), 1, 1, length(theta)),reshape(sin(theta),1,1,[]));
    exdim = repmat(qvsintheta, 1, 1, 1, length(phi)); %add extra dimension over phi
    B = bsxfun(@times,exdim,reshape(cos(phi),1,1,1,[]));
    
    w_wave = bsxfun(@minus,repmat((osc.eloss/h2ev)',1,2^(decdigs-1)+1,length(theta),length(phi)),B.*sind(alpha));
    Qv_per = Q.^2.*v_per^2;
    bottom = bsxfun(@plus,w_wave.^2,Qv_per);
    
    exdimdep = repmat(Q, 1, 1, 1, length(phi),length(r)); %add extra dimension over r
    exprQ = bsxfun(@times,reshape((-1)*abs(r),1,1,1,1,[]),exdimdep).*cosd(alpha);
    top_in = bsxfun(@times,bsxfun(@times,qsintheta,cos(qzrcos)),exp(exprQ));
    topbot = bsxfun(@rdivide,top_in,repmat(bottom, 1, 1, 1, 1, length(r)));
    ind0 = bsxfun(@eq,depth,0);
    ind = bsxfun(@gt,depth,0);
    
    %% Surface
    x_in = zeros(length(osc.eloss),length(depth));
    
    for i = 1:length(theta)
        osc.qtran = Q(:,:,i)/a0;
        Im(:,:,i) = eps_sum_allwq(osc,'surface');
    end
    
    %================= inside ===================
    romall_in = bsxfun(@times,Im,topbot);
    romall_in(isnan(romall_in)) = 0;
    romall_in(:,:,:,:,ind0) = romall_in(:,:,:,:,ind0).*0.5;
    romall_in(:,:,:,:,ind) = 0;
    res_in = squeeze(trapz(theta,trapz(phi,romall_in,4),3));

    
    %================= outside ==================
    top_out = bsxfun(@times,qsintheta,exp(exprQ));
    exdimdep = repmat(w_wave, 1, 1, 1, 1,length(r)); %add extra dimension over r
    cosw = bsxfun(@times,exdimdep,reshape(r,1,1,1,1,[]));
    add = bsxfun(@minus,2.*cos(cosw./sqrt(2*E0/h2ev)),exp(exprQ));
    
    romall_out = bsxfun(@times,bsxfun(@times,Im,bsxfun(@rdivide,top_out,bottom)),add);
    romall_out(isnan(romall_out))=0;
    romall_out(:,:,:,:,ind0) = romall_out(:,:,:,:,ind0).*0.5;
    ind = bsxfun(@lt,depth,0);
    romall_out(:,:,:,:,ind) = 0;
    res_out = squeeze(trapz(theta,trapz(phi,romall_out,4),3));
    
    for i=1:length(osc.eloss)
        for j = 1:length(depth)
            x_in(i,j) = 4*cosd(alpha)/(pi^3)/h2ev/a0 * (trapz(q(i,:),res_in(i,:,j)) + trapz(q(i,:),res_out(i,:,j)));
        end
    end
    
    siimfp = trapz(osc.eloss,x_in);
   
    dsep = x_in; %./trapz(osc.eloss,x_in); % only surface component 
end


