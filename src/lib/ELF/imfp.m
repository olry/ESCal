function l_in = imfp(osc,energy,material,Penn,elf,varargin)

    if nargin < 3
        material = 'metal';
    end
    if nargin < 4
        Penn = false;
    end
    
    if energy<500
        osc.eloss = fliplr(energy:-.5:eps);
    else
        osc.eloss = [fliplr(500:-.5:eps) fliplr(energy:-10:510)];
    end
    if Penn
        disp('Using the Penn algorithm ...');
        if nargin < 5
            elf = eps_sum(osc);
            l_in = iimfp_penn(energy,osc.eloss,elf);
        else
            l_in = iimfp_penn(energy,elf(:,1),elf(:,2));
        end            
    else
        w_in = ndiimfp(osc,energy,11,false);
        if strcmp(material,'metal')
            int_emesh = eps:0.1:energy - osc.Ef;
        elseif strcmp(material,'insulator')
            int_emesh = eps:0.1:energy - (osc.egap + osc.vb);
        end
        xin_new = interp1(osc.eloss,w_in,int_emesh);
        l_in = 1/trapz(int_emesh,xin_new);
    end  
end