function [energy,elf] = mopt(osc,path2xray,merge)

    if nargin<3
        merge = false;
    end

    f = load([path2xray osc.formula.symbols{1} '.']);
    num = length(osc.formula.symbols);
    optData = zeros(length(osc.formula.symbols),length(f),3);
    f1sum = 0;
    f2sum = 0;

    for i = 1:num
        optData(i,:,:) = load([path2xray osc.formula.symbols{i} '.']);
        f1sum = f1sum + optData(i,:,2)*osc.formula.weights(i);
        f2sum = f2sum + optData(i,:,3)*osc.formula.weights(i);
    end
    
    f1sum = f1sum/sum(osc.formula.weights);
    f2sum = f2sum/sum(osc.formula.weights);
          
    lambda = hc./(optData(i,:,1)./1000);
    n = 1 - osc.na*r0*1e10*(lambda.^2).*(f1sum)/2/pi;
    k = -osc.na*r0*1e10*(lambda.^2).*(f2sum)/2/pi;
    
    eps1 = n.^2 - k.^2;
    eps2 = 2.*n.*k;
    
    if merge
        osc.eloss = eps:0.01:100;
        ind = bsxfun(@gt,optData(i,:,1),100); 
        elf = [eps_sum(osc,true)' -eps2(ind)./(eps1(ind).^2 + eps2(ind).^2)];
        energy = [osc.eloss optData(i,ind,1)];
    else
        elf = -eps2./(eps1.^2 + eps2.^2);
        energy = optData(i,:,1);
    end
end

