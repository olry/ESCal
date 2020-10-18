function osc_scaled = scaling(osc,xraypath,varargin)

if nargin<2
    xraypath = '';
end

    o = convert2au(osc);
    if xraypath
        [eloss,elf_henke] = mopt(osc,xraypath,false);
        ind = bsxfun(@gt,eloss,100);
        bsum_henke = 1/(2*pi^2)*trapz(eloss(ind)/h2ev,bsxfun(@times,eloss(ind)/h2ev,elf_henke(ind)));
    else
        bsum_henke = 0;
    end
    switch o.model
        case 'Drude'
            o.A = o.A/sum(o.A)*4*pi*(o.ne*a0^3 - bsum_henke);
%             eps1 = eps_real_Drude(o);
%             w_ScalingFactor = sqrt((eps1(1) - o.beps)/(o.n_refrac^2 - o.beps));
%             o.Om = o.Om * w_ScalingFactor;
        case 'DrudeLindhard'
            o.A = o.A/sum(o.A)*(1-1/o.n_refrac^2);
            BetheSum = sum((pi/2)*o.A.*o.Om.^2);
            BetheValue = 2*pi*pi*(o.ne*a0^3 - bsum_henke); %o.ne*a0^3;
            w_ScalingFactor = sqrt(BetheSum / BetheValue);
            o.Om = o.Om / w_ScalingFactor;
        case 'Mermin'
            o.A = o.A/sum(o.A)*(1-1/o.n_refrac^2);
            BetheSum = sum((pi/2)*o.A.*o.Om.^2);
            BetheValue = 2*pi*pi*(o.ne*a0^3 - bsum_henke); %o.ne*a0^3;
            w_ScalingFactor = sqrt(BetheSum / BetheValue);
            o.Om = o.Om / w_ScalingFactor;
        case 'MerminLL'
            o.A = o.A/sum(o.A);
            BetheSum = sum((pi/2)*o.A.*o.Om.^2);
            BetheValue = 2*pi*pi*(o.ne*a0^3 - bsum_henke); %o.ne*a0^3;
            w_ScalingFactor = sqrt(BetheSum / BetheValue);
            o.Om = o.Om / w_ScalingFactor;
%             u_step = 0.5/h2ev;
%             enlargeU = false;
% 			decreaseU = false;
%             eps1 = eps_real_MLL(o);
%             while abs(eps1 - o.n_refrac^2) > 0.001
%                 if (eps1 - o.n_refrac^2) > 0.001
%                     o.u = o.u + u_step;
% 					enlargeU = true;
%                 elseif (eps1 - o.n_refrac^2) < 0.001
%                     o.u = o.u - u_step;
%                     decreaseU = true;
%                 end
%                 if decreaseU && enlargeU
%                     decreaseU = false;
% 					enlargeU = true;
% 					u_step = u_step / 5;
%                 end
%                 eps1 = eps_real_MLL(o);
%             end
            
    end
    osc_scaled = convert2ru(o);


    function eps = eps_real_Drude(o)
        w = o.eloss;
        q = o.qtran;
        eps_re = o.beps;
        for j=1:length(o.A)
            epsDrud_re = Drude(q,w,o.Om(j),o.G(j),o.alpha,o.Ef);
            eps_re = eps_re - o.A(j)*epsDrud_re;
        end
        eps = eps_re;
    end

    function eps_1 = eps_real_MLL(o)
        w = 0.0001;
        q = 0.01;
        eps = 0;
        for j=1:length(o.A)
            epsMerm = Mermin_LL(q,w,o.G(j),o.Om(j),o.u);
            eps = eps + o.A(j)*(complex(1,0)./epsMerm);
        end
        eps_tot = complex(1,0)./eps;
        eps_1 = real(eps_tot);
    end

end
