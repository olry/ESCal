function fit_result = fit_xps_spectrum(data,osc,xraypath)

eval_num = 0;

%% constraints
osc_min.A = ones(size(osc.A))*1e-10;
osc_min.G = ones(size(osc.G))*0.02; 
osc_min.Om = ones(size(osc.Om));

switch osc.model
    case 'Drude'
        coef = 1e3;
        osc_max.G = ones(size(osc.G))*100;
    case 'DrudeLindhard'
        coef = 1;
        osc_max.G = ones(size(osc.G))*100;
    case 'Mermin'
        coef = 1;
        osc_max.G = ones(size(osc.G))*30;
    case 'MerminLL'
        coef = 1;
        osc_max.G = ones(size(osc.G))*30;
        osc_min.Om = zeros(size(osc.Om));
        osc_min.u = 0;
        osc_max.u = 20;
end
       
osc_max.A = ones(size(osc.A))*coef;
osc_max.Om = ones(size(osc.Om))*100;

lb = structToVec(osc_min);
ub = structToVec(osc_max);

%% NLopt
% {
opt.algorithm = NLOPT_LN_COBYLA;
opt.lower_bounds = lb;
opt.upper_bounds = ub;
opt.maxeval = 100;
opt.min_objective = @fit_func_nlopt;
opt.fc = { (@(x) constraint(x)) };
opt.fc_tol = 1e-8; 
opt.xtol_rel = 1e-10;
[x_res] = nlopt_optimize(opt, structToVec(osc));
an = vecToStruct(x_res);
%}

disp(an);
fit_result = an;

%% Plotting
figure;
plot(data.x_exp,data.y_exp,'DisplayName','Experiment','Marker','o','LineWidth',1)
hold on
plot(data.x_exp,fit_func(x_res,data.x_exp),'DisplayName','Summary signal','LineWidth',2)
legend

xlabel('Kinetic energy, eV')
ylabel('Intensity, rel.un.')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% Sum-rules

an_au = convert2au(an);
[eloss,elf] = mopt(an,xraypath,true);

bsum = 1/(2*pi^2)*trapz(eloss/h2ev,bsxfun(@times,eloss/h2ev,elf));
psum = 2/pi*trapz(eloss(2:end),bsxfun(@rdivide,elf(2:end),eloss(2:end))) + 1/an.n_refrac^2;
fsum = 1/(2*pi^2*(an.na*a0^3))*trapz(eloss/h2ev,bsxfun(@times,eloss/h2ev,elf));

disp(['P-sum rule: ',num2str(psum)]);
disp(['Bethe sum rule: ',num2str(bsum), ' electron density = ',num2str(osc.ne*a0^3), '(a.u.^-3)']);
disp(['Sum of A: ',num2str(sum(an_au.A)),'(1-1/n2) = ', num2str(abs(1-(1/osc.n_refrac)^2)), ' 4piN = ',num2str(4*pi*osc.ne*a0^3), '(a.u.^-3)']);
disp(['F-sum rule:',num2str(fsum), ' EAN = ', num2str(osc.Z)]);

%%
function v = structToVec(s)
    v = [s.A, s.G, s.Om];
end

function s = vecToStruct(v)
    s = osc;
    nA = length(osc.A);

    s.A = v(1:nA);
    s.G = v((1:nA)+nA);
    s.Om = v((1:nA)+nA*2);
end

function y = fit_func(os,xdata)

    o = vecToStruct(os);

    %% x_in        
    diimfp = ndiimfp(o,data.E0,10,true,false);
    x_in_b = interp1(o.eloss,diimfp,data.mesh_eloss);
    x_in_b(isnan(x_in_b)) = eps;

    %% spectrum bulk
    data.Mat.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_b'));
    sigma_lds = 0.08; 
    alpha_lds = 0.04;
    MS.C = {'1S1/2'};
    Layers = Layer(data.Mat);
    Q = PESMultiLayer(Layers,MS);
    Q.N_in = 15;
    Q.theta0 = data.theta0;
    Q.theta = data.theta;
    Q.phi = 180 - data.phi;
    Q.energy_mesh_full = xdata';
    Q.sigma_gauss = data.sigma_C; 
    Q.alpha_LDS = alpha_lds;
    Q.sigma_LDS = sigma_lds;
    Q.Calculate;
    Q.CalculateEnergyDistribution(data.theta,180-data.phi);
    delta = 3.3;
    delta_2 = 2;
    CO_peak = interp1(Q.energy_mesh_full,Q.EnergyDistribution,Q.energy_mesh_full+delta,'spline','extrap');
    CO_2_peak = interp1(Q.energy_mesh_full,Q.EnergyDistribution,Q.energy_mesh_full+delta_2,'spline','extrap');
    res = Q.EnergyDistribution/max(Q.EnergyDistribution) + CO_peak/max(CO_peak)*0.015 + CO_2_peak/max(CO_2_peak)*0.015;
    res = (res+0.005);
    y = interp1(data.x_exp,res/max(res),data.x_exp-2.68);  
    y(isnan(y))=eps;
end


function [val, gradient] = fit_func_nlopt(os)
    eval_num = eval_num + 1;
    o = vecToStruct(os);

    %% x_in        
    diimfp = ndiimfp(o,data.E0,10,true,false);
    x_in_b = interp1(o.eloss,diimfp,data.mesh_eloss);
    x_in_b(isnan(x_in_b)) = eps;

    %% spectrum bulk
    data.Mat.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_b'));
    sigma_lds = 0.08; 
    alpha_lds = 0.04;
    MS.C = {'1S1/2'};
    Layers = Layer(data.Mat);
    Q = PESMultiLayer(Layers,MS);
    Q.N_in = 15;
    Q.theta0 = data.theta0;
    Q.theta = data.theta;
    Q.phi = 180 - data.phi;
    Q.energy_mesh_full = data.x_exp';
    Q.sigma_gauss = data.sigma_C; 
    Q.alpha_LDS = alpha_lds;
    Q.sigma_LDS = sigma_lds;
    Q.Calculate;
    Q.CalculateEnergyDistribution(data.theta,180-data.phi);
    delta = 3.3;
    delta_2 = 2;
    CO_peak = interp1(Q.energy_mesh_full,Q.EnergyDistribution,Q.energy_mesh_full+delta,'spline','extrap');
    CO_2_peak = interp1(Q.energy_mesh_full,Q.EnergyDistribution,Q.energy_mesh_full+delta_2,'spline','extrap');
    res = Q.EnergyDistribution/max(Q.EnergyDistribution) + CO_peak/max(CO_peak)*0.015 + CO_2_peak/max(CO_2_peak)*0.015;
%     res = (res+0.005);
    y = interp1(data.x_exp,res/max(res),data.x_exp-2.68);  
    y(isnan(y))=eps;
      
    val = sum((data.y_exp - y).^2);
    if mod(eval_num,1) == 0
        disp(['Number of function evaluations:', num2str(eval_num), ' chisq:', num2str(val)]);
    end
    if (nargout > 1)
        gradient = [0, 0.5 / val];
    end
end

function [val, gradient] = constraint(x)
    o = vecToStruct(x);
    o_au = convert2au(o);

    [eloss,elf_henke] = mopt(o,xraypath);
    ind = bsxfun(@gt,eloss,100);
    bsum_henke = 1/(2*pi^2)*trapz(eloss(ind)/h2ev,bsxfun(@times,eloss(ind)/h2ev,elf_henke(ind)));
    cf = 4*pi*(o.ne*a0^3 - bsum_henke) / sum(o_au.A);         
	val = abs(cf-1);

    if (nargout > 1)
        gradient = [0, 0.5 / val];
    end   
end

end
