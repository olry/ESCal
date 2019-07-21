function [diimfp,dsep,siimfp] = ndiimfp_Li(osc,E0,depth,alpha,decdigs,varargin)

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
       
    qmin = log(sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-osc.eloss/h2ev)));
    qmax = log(sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-osc.eloss/h2ev)));
    q = zeros(length(osc.eloss),2^(decdigs-1)+1);
    
    for i = 1:2^(decdigs-1)+1
        q(:,i) = qmin + (i-1)*(qmax-qmin)/2.^(decdigs-1);
    end
    
    if strcmp(osc.model,'Mermin')
        q(q==0) = 0.01;
    end
    
    osc.qtran = exp(q)/a0;
    
    %% Clear bulk
    x_in_b = zeros(size(osc.eloss));
    
    ELF = eps_sum_allwq(osc,'bulk');
    ELF(isnan(ELF)) = 0;

    x_in_b(1) = 0;
    for i = 2:length(osc.eloss)
        x_in_b(i) = 1/pi/(E0/h2ev) * trapz(q(i,:),ELF(i,:))/h2ev/a0;
    end

    %% Surface
    
    qmin = sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-osc.eloss/h2ev));
    qmax = sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-osc.eloss/h2ev));

    for i = 1:2^(decdigs-1)+1
        q(:,i) = qmin + (i-1)*(qmax-qmin)/2.^(decdigs-1);
    end
    
    if strcmp(osc.model,'Mermin')
        q(0.001<q==0) = 0.01;
    end
    osc.qtran = q/a0;
        
    theta = 0:pi/2/10:pi/2;
    phi = 0:2*pi/10:2*pi;
    
    Im = zeros(length(osc.eloss),2^(decdigs-1)+1,length(theta));
    x_in = zeros(length(osc.eloss),length(depth));
    
    Q = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(sin(theta),1,1,[]));
    indQ = bsxfun(@eq,Q,0);
    Q(indQ) = 0.0001;
       
    qz = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(cos(theta),1,1,[]));
    v_per = cosd(alpha).*sqrt(2*E0/h2ev);
    r = depth./a0./cosd(alpha);
    %exdimr = repmat(qz, 1, 1, 1,length(phi)); %add extra dimension over phi
%     qzrcos = bsxfun(@times,exdimr,cosd(alpha));
    qzrcos = qz.*cosd(alpha);
    
    qsintheta = bsxfun(@times,repmat(q, 1, 1, length(theta)),reshape(sin(theta).^2,1,1,[]));
    
    qvsintheta = bsxfun(@times,repmat(q.*sqrt(2*E0/h2ev), 1, 1, length(theta)),reshape(sin(theta),1,1,[]));
    exdim = repmat(qvsintheta, 1, 1, 1, length(phi)); %add extra dimension over phi
    B = bsxfun(@times,exdim,reshape(cos(phi),1,1,1,[]));
    
    w_wave = bsxfun(@minus,repmat((osc.eloss/h2ev)',1,2^(decdigs-1)+1,length(theta),length(phi)),B.*sind(alpha));
    Qv_per = Q.^2.*v_per^2;
    bottom = bsxfun(@plus,w_wave.^2,Qv_per);
              
    for i = 1:length(theta)
        osc.qtran = Q(:,:,i)/a0;
        Im(:,:,i) = eps_sum_allwq(osc,'surface');
    end
    
    for dep = 1:length(r)
        top_in = (qsintheta.*cos(qzrcos.*r(dep))).*exp(Q.*cosd(alpha).*(-abs(r(dep))));   
        topbot = bsxfun(@rdivide,top_in,bottom);
        %================= inside ===================
        romall_in = bsxfun(@times,Im,topbot);
        romall_in(isnan(romall_in)) = 0;
        romall_in = romall_in.*stepfunction(r(dep));
        res_in = squeeze(trapz(theta,trapz(phi,romall_in,4),3));

        %================= outside ==================
        top_out = bsxfun(@times,qsintheta,exp(Q.*cosd(alpha).*(-abs(r(dep)))));
        cosw = w_wave.*r(dep);
        add = 2.*cos(cosw./sqrt(2*E0/h2ev))-exp(Q.*cosd(alpha).*(-abs(r(dep))));

        romall_out = bsxfun(@times,Im,top_out./bottom).*add;
        romall_out(isnan(romall_out))=0;
        romall_out = romall_out.*stepfunction(-r(dep));
        res_out = squeeze(trapz(theta,trapz(phi,romall_out,4),3));

        for i=1:length(osc.eloss)
            x_in(i,dep) = 4*cosd(alpha)/(pi^3)/h2ev/a0 * (trapz(q(i,:),res_in(i,:)) + trapz(q(i,:),res_out(i,:)));
        end   
    end
    
    %% siimfp/biimfp
    siimfp = trapz(osc.eloss,x_in); % + x_b + x_in_b,1);
%     biimfp = trapz(osc.eloss,x_in_b+x_b,1);
%     biimfp(biimfp<0)=0;
%     iimfp = trapz(osc.eloss,x_in_b+x_b+x_in,1);
%     figure;
%     plot(depth,siimfp,depth,biimfp,depth,siimfp+biimfp)
    
    %% Plot
%     figure;
%     xlim([0 100])
%     hold on
%     box on
%     plot(osc.eloss,x_b(:,1))
%     plot(osc.eloss,x_b(:,1) + x_in_b(:,1))
%     plot(osc.eloss,x_in(:,1))
%     plot(osc.eloss,x_b(:,1) + x_in_b(:,1) + x_in(:,1))
%     legend('Clear bulk','Reduced bulk','Surface','DIIMFP');
%     
%     Y = ['siimfp = ',num2str(siimfp)];
%     disp(Y);
%     Y = ['biimfp = ',num2str(biimfp)];
%     disp(Y);
    
    %dsep = (x_in+x_b+x_in_b)./trapz(osc.eloss,x_in+x_b+x_in_b);
    %dsep = x_in + x_b + x_in_b; %./trapz(osc.eloss,x_in); 
    dsep = x_in; %./trapz(osc.eloss,x_in); % only surface component
    diimfp = x_in_b; %./trapz(osc.eloss,x_b);   % clear bulk  
end

%% Heaviside function
    function x = stepfunction(depth)
        if depth > 0
            x = 1;
        elseif depth < 0
            x = 0;
        elseif depth==0
            x = 0.5;
        end
    end

end










