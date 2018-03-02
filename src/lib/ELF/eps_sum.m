function ELF = eps_sum(osc)

%%
%{
   Calculates the imaginary part of the inverse dielectric function Im[-1/eps]
   for the energy w and momentum transfer q on the basis of four ELF models.
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

if strcmp(osc.model,'Drude') || strcmp(osc.model,'DrudeLindhard')
    osc = convert2au(osc); %converts to atomic units, this is necessary for Drude model
end

elec_density = sum(osc.A)/4/pi;
Y = ['Electron density is ',num2str(elec_density),' (in a.u.)'];
%disp(Y);

w = osc.eloss;
q = osc.qtran;

ELF = zeros(numel(w),numel(q));

if strcmp( osc.model,'Drude')
    eps_re = osc.beps*ones(numel(w),numel(q));
    eps_im = zeros(numel(w),numel(q));
    for j=1:length(osc.A)      
        
        
        if isfield(osc,'numion') && j>length(osc.A) - osc.numion
            [epsDrud_re, epsDrud_im] = Drude(q,w,osc.Om(j),osc.G(j),osc.alpha,osc.Ef,true);
        else
            [epsDrud_re, epsDrud_im] = Drude(q,w,osc.Om(j),osc.G(j),osc.alpha,osc.Ef,false);
        end
        eps_re = eps_re - osc.A(j)*epsDrud_re;
        eps_im(w>osc.egap,:) = eps_im(w>osc.egap,:) + osc.A(j)*epsDrud_im(w>osc.egap,:);
        %plot(w,imag(-1./complex(osc.beps-osc.A(j)*epsDrud_re,osc.A(j)*epsDrud_im(w>osc.egap,:))));
        %plot(w,imag(-1./complex(eps_re,eps_im)));
    end
    eps = complex(eps_re,eps_im);
    ELF = imag(-1./eps);
elseif strcmp( osc.model,'DrudeLindhard')
    
    eps_re = ones(numel(w),numel(q));
    eps_im = zeros(numel(w),numel(q));
    for j=1:length(osc.A)      
        
        
        if isfield(osc,'numion') && j>length(osc.A) - osc.numion
            [epsDrud_re, epsDrud_im] = DrudeLindhard(q,w,osc.Om(j),osc.G(j),osc.alpha,osc.Ef,true);
        else
            [epsDrud_re, epsDrud_im] = DrudeLindhard(q,w,osc.Om(j),osc.G(j),osc.alpha,osc.Ef,false);
        end
        eps_re = eps_re + osc.A(j)*epsDrud_re;
        eps_im(w>osc.egap,:) = eps_im(w>osc.egap,:) + osc.A(j)*epsDrud_im(w>osc.egap,:);
        %plot(w,imag(-1./complex(osc.beps-osc.A(j)*epsDrud_re,osc.A(j)*epsDrud_im(w>osc.egap,:))));
        %plot(w,imag(-1./complex(eps_re,eps_im)));
        %plot(w,eps_im);
    end
    eps = complex(eps_re,eps_im);
    %ELF = imag(-1./eps);
    ELF = eps_im;
else

    for k=1:length(q)
        for i=1:length(w)
            eps1 = complex(0.0,0.0);
            for j=1:length(osc.A)      
                switch osc.model
                    case 'Lindhard'
                        w_q = osc.Om(j) + 0.5*osc.alpha * q(k)*q(k);
                        epsLind = Lindhard(q(k),w(i),osc.G(j),w_q);
                        eps1 = eps1 + osc.A(j)*(complex(1,0)/epsLind);
                    case 'Mermin'
                        w_q = osc.Om(j) + 0.5*osc.alpha * q(k)*q(k);
                        epsMerm = Mermin(q(k),w(i),osc.G(j),w_q);
                        eps1 = eps1 + osc.A(j)*(complex(1,0)/epsMerm);
                    case 'Mermin_LL'
                        w_q = osc.Om(j) + 0.5*osc.alpha * q(k)*q(k);
                        epsMerm = Mermin_LL(q(k),w(i),osc.G(j),w_q,osc.u);
                        eps1 = eps1 + osc.A(j)*(complex(1,0)/epsMerm);
                    otherwise
                        error('Choose the correct model name: Drude,Lindhard,Mermin,Mermin_LL');
                end                
            end 
            eps = complex(1,0)/eps1;

            ELF(i,k) = imag(-1/eps);
        end
    end 
end
end