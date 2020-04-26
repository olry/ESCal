function [energy,elf] = mopt(osc)

    f = load(['/Users/olgaridzel/Bruce/PHYSDAT/opt/xray/' osc.formula.symbols(1) '.']);
    optData = zeros(length(osc.formula.symbols),length(f),3);
    f1sum = 0;
    f2sum = 0;

    for i = 1:length(osc.formula.symbols)
        optData(i,:,:) = load(['/Users/olgaridzel/Bruce/PHYSDAT/opt/xray/' osc.formula.symbols(i) '.']);
        f1sum = f1sum + optData(i,:,2)*osc.formula.weights(i);
        f2sum = f2sum + optData(i,:,3)*osc.formula.weights(i);
    end
    
    f1sum = f1sum/sum(osc.formula.weights);
    f2sum = f2sum/sum(osc.formula.weights);
    
    osc.eloss = eps:0.01:100;
    ind = bsxfun(@gt,optData(i,:,1),100);    
    energy = [osc.eloss optData(i,ind,1)];
       
    lambda = hc./(optData(i,:,1)./1000);
    n = 1 - osc.na*r0*1e10*(lambda.^2).*(f1sum)/2/pi;
    k = -osc.na*r0*1e10*(lambda.^2).*(f2sum)/2/pi;
    
    eps1 = n.^2 - k.^2;
    eps2 = 2.*n.*k;
    
    elf = [eps_sum(osc)' -eps2(ind)./(eps1(ind).^2 + eps2(ind).^2)];
end