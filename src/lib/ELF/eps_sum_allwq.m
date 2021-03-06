function ELF = eps_sum_allwq(osc,interface,gapEffect,varargin)

if nargin < 3
    gapEffect = false;
end

%%
%{
   Calculates the imaginary part of the inverse dielectric function Im[-1/eps]
   for the energy array w and momentum transfer array q on the basis of four ELF models
   for the consequent integration over all momentum transfers.
   Input parameter:
   osc is the field of oscilattors data including the follows:
         1. osc.A      - amplitudes
         2. osc.G      - damping coefficients gamma
         3. osc.Om     - plasmon (resonance) frequency omega
         4. osc.alpha  - alpha is a constant between 1 (metals) and 0 (insulators)
         5. osc.u      - a quantity related to the band
                        gap (is used only for Mermin_LL)
         6. osc.model  - type of model: Drude, Lindhard, Mermin, Mermin_LL
         7. osc.eloss  - energy loss array (e.g. 0:0.5:100 eV)
         8. osc.qtran  - momentum transfer array (e.g. 0:0.01:10 A^-1)
                        In the case of Mermin and Mermin_LL models the array of q should always
                        start from 0.01 (instead of 0).
         9. osc.Ef     - the Fermi energy
        10. osc.beps   - the background dielectric constant due to the polarizability of the core electrons
                        osc.beps = 1 assuming all electrons are free
        11. osc.egap   - the band gap energy (for insulators)
    
    For more details see Vos, M., and Grande, P. L. (2017) Extracting the dielectric function from high-energy REELS measurements. Surf. Interface Anal., 49: 809–821. doi: 10.1002/sia.6227
%}
%%

osc = convert2au(osc); %converts to atomic units

w = osc.eloss;
q = osc.qtran;

if strcmp( osc.model,'Drude')
    
    eps_re = osc.beps*ones(size(q));
    eps_im = zeros(size(q));
    for j=1:length(osc.A)
        [epsDrud_re, epsDrud_im] = Drude(q,w,osc.Om(j),osc.G(j),osc.alpha,osc.Ef);
        eps_re = eps_re - osc.A(j)*epsDrud_re;
        if gapEffect
            ind = bsxfun(@gt,w,osc.egap);
            eps_im(ind,:) = eps_im(ind,:) + osc.A(j)*epsDrud_im(ind,:);
        else
            eps_im = eps_im + osc.A(j)*epsDrud_im;
        end
    end
    epsilon = complex(eps_re,eps_im);
    if strcmp(interface,'bulk')
%         ELF = imag(-1./epsilon);
        ELF = eps_im ./ (eps_re.^2 + eps_im.^2);
    elseif strcmp(interface,'surface')
%         den = (eps_re.^2 + eps_re - eps_im.^2).^2 + (2*eps_re.*eps_im + eps_im).^2;
%         enu = -eps_im.*(2*eps_re + 1).*((eps_re - 1).^2 - eps_im.^2);
%         enu = enu + 2*eps_im.*(eps_re - 1).*(eps_re.*(eps_re + 1) - eps_im.^2);
%         ELF = enu./den;
        ELF = imag(-1./(epsilon + 1));
    end
elseif strcmp( osc.model,'DrudeLindhard')
    sumoneovereps = complex(1,0);   
    for j=1:length(osc.A)
        [epsDrud_re, epsDrud_im] = DrudeLindhard(q,w,osc.Om(j),osc.G(j),osc.alpha,osc.Ef);
        oneoverepsDL = complex(epsDrud_re, -epsDrud_im);
        sumoneovereps = sumoneovereps + osc.A(j) * (oneoverepsDL - complex(1, 0));
    end
    epsilon = complex(1, 0) ./ sumoneovereps;
    if strcmp(interface,'bulk')
        ELF = imag(-1./epsilon);
    elseif strcmp(interface,'surface')
        eps_re = real(epsilon);
        eps_im = imag(epsilon);
        den = (eps_re^2 + eps_re - eps_im^2)^2 + (2*eps_re*eps_im + eps_im)^2;
        enu = -eps_im*(2*eps_re + 1)*((eps_re - 1)^2 - eps_im^2);
        enu = enu + 2*eps_im*(eps_re - 1)*(eps_re*(eps_re + 1) - eps_im^2);
        ELF = enu/den;
    end
    if gapEffect
        ELF(w<osc.egap,:) = 0;
    end
  
%     eps_re = zeros(size(q));
%     eps_im = zeros(size(q));
%     for j=1:length(osc.A)
%         [epsDrud_re, epsDrud_im] = DrudeLindhard(q,w,osc.Om(j),osc.G(j),osc.alpha,osc.Ef);
%         eps_re = eps_re + osc.A(j)*epsDrud_re;
%         ind = bsxfun(@gt,w,osc.egap);
%         eps_im(ind,:) = eps_im(ind,:) + osc.A(j)*epsDrud_im(ind,:);
%     end
%     if strcmp(interface,'bulk')
%         if xraypath
%             ELF = get_henke_data(eps_im);
%         else
%             ELF = eps_im;
%         end
%     elseif strcmp(interface,'surface')
%         error('The DrudeLindhard model cannot be used for surface calculations');
%     end
elseif strcmp( osc.model,'Mermin')
    eps1 = zeros(size(q));
    for j=1:length(osc.A)
        epsMerm = Mermin(q,w,osc.G(j),osc.Om(j));
        eps1 = eps1 + osc.A(j)*(complex(1,0)./epsMerm);
    end
    epsilon = complex(1,0)./eps1;
    if strcmp(interface,'bulk')
        ELF = imag(-1./epsilon);
    elseif strcmp(interface,'surface')
        eps_re = real(epsilon);
        eps_im = imag(epsilon);
        den = (eps_re^2 + eps_re - eps_im^2)^2 + (2*eps_re*eps_im + eps_im)^2;
        enu = -eps_im*(2*eps_re + 1)*((eps_re - 1)^2 - eps_im^2);
        enu = enu + 2*eps_im*(eps_re - 1)*(eps_re*(eps_re + 1) - eps_im^2);
        ELF = enu/den;
    end
    if gapEffect
        ELF(w<osc.egap,:) = 0;
    end
elseif strcmp( osc.model,'MerminLL')
    eps1 = zeros(size(q));
    for j=1:length(osc.A)
        epsMerm = Mermin_LL(q,w,osc.G(j),osc.Om(j),osc.u);
        eps1 = eps1 + osc.A(j)*(complex(1,0)./epsMerm);
    end
    epsilon = complex(1,0)./eps1;
    ELF = imag(-1./epsilon);
else
    error('Choose the correct model name: Drude,DrudeLindhard,Mermin,MerminLL');
    
end

end