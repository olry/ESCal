function fit_result = optfit_spectrum(data,osc,xraypath)

eval_num = 0;

%% constraints
osc_min.A = ones(size(osc.A))*1e-10;
osc_min.G = ones(size(osc.G))*0.02; 
osc_min.Om = ones(size(osc.Om));
osc_min.egap = 1.7;
osc_min.H = eps;
osc_min.D = eps;

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
osc_max.egap = 3.2;
osc_max.H = 0.06;
osc_max.D = 0.06;

lb = structToVec(osc_min);
ub = structToVec(osc_max);

%% local minimum search
%{
options = optimoptions('lsqcurvefit','PlotFcn',@optimplotx,'Display','iter-detailed','UseParallel',true);
% options.StepTolerance = 1e-12;
% options.MaxFunctionEvaluations = 5000;
[x_res] = lsqcurvefit(@fit_func, structToVec(osc),data.x_exp,data.y_exp,lb,ub,options);
an = vecToStruct(x_res);
an = scaling(an,xraypath);
%}

%% NLopt
% {
opt.algorithm = NLOPT_LN_COBYLA;
opt.lower_bounds = lb;
opt.upper_bounds = ub;
opt.maxeval = 2000;
opt.min_objective = @fit_func_nlopt;
opt.fc = { (@(x) aconstraint(x)) };
opt.fc_tol = 1e-8; 
opt.xtol_rel = 1e-10;
[x_res] = nlopt_optimize(opt, structToVec(osc));
an = vecToStruct(x_res);
%}

disp(an);
fit_result = an;

%% Plotting

% plot(data.E0 - data.x_exp,data.y_exp,'DisplayName','Experiment','Marker','o','LineWidth',1)
% hold on
% plot(data.E0 - data.mesh_eloss,data.Gauss_H*data.int_H,'DisplayName','H signal','LineWidth',2)
% plot(data.E0 - data.mesh_eloss,data.Gauss_D*data.int_D,'DisplayName','D signal','LineWidth',2)
% plot(data.E0 - data.x_exp,fit_func(x_res,data.x_exp),'DisplayName','Summary signal','LineWidth',2)

figure;
plot(data.E0 - data.x_exp,data.y_exp,'DisplayName','Experiment','Marker','o','LineWidth',1)
hold on
plot(data.E0 - data.x_exp,fit_func(x_res,data.x_exp),'DisplayName','Summary signal','LineWidth',2)
legend

xlabel('Kinetic energy, eV')
ylabel('Intensity, rel.un.')
set(findall(gcf,'-property','FontSize'),'FontSize',25)

% ylim([0 0.05])
% xlim([-3 40])

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% txt = 'chsurfs';
% print(h,txt,'-dpdf','-r0')
% txt = 'chsurfs.fig';
% savefig(txt)

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
    v = [s.A, s.G, s.Om, s.egap, s.H, s.D];
%     v = [s.A, s.G, s.Om];
end

function s = vecToStruct(v)
    s = osc;
    nA = length(osc.A);

    s.A = v(1:nA);
    s.G = v((1:nA)+nA);
    s.Om = v((1:nA)+nA*2);
    s.egap = v(end-2);
    s.H = v(end-1);
    s.D = v(end);
end

function [x_in_b, x_in_s, int_over_depth_sigma_surf] = crosssection(o)
    [diimfp,dsep,sigma] = ndiimfp_Li(o,data.E0,data.depth,data.theta,8,xraypath);
    r = data.depth./cosd(data.theta);

    int_over_depth_dsep = trapz(r,dsep,2);
    int_over_depth_sigma_surf = trapz(r,sigma); % SEP

    xs_new = interp1(o.eloss,int_over_depth_dsep,data.mesh_eloss);
    xs_new(isnan(xs_new)) = 0;
    x_in_s = xs_new./trapz(data.mesh_eloss,xs_new); % dsep

    xb_new = interp1(o.eloss,diimfp,data.mesh_eloss);
    xb_new(isnan(xb_new)) = 0;
    x_in_b = xb_new./trapz(data.mesh_eloss,xb_new); % diimfp
end

function y = fit_func(os,xdata)

    o = vecToStruct(os);
%     o = scaling(o,xraypath);

    %% x_in        
%     [x_in_b, x_in_s, int_over_depth_sigma_surf] = crosssection(o);
    dsep = dsep_Tung_2(o, data.E0, 12, false, 60);
    diimfp = ndiimfp(o, data.E0, 12, false);

    sep = trapz(o.eloss,dsep); %*h2ev;

    xs_new = interp1(o.eloss,dsep,data.mesh_eloss);
    xs_new(isnan(xs_new)) = 0;
    x_in_s = xs_new./trapz(data.mesh_eloss,xs_new); % dsep

    xb_new = interp1(o.eloss,diimfp,data.mesh_eloss);
    xb_new(isnan(xb_new)) = 0;
    x_in_b = xb_new./trapz(data.mesh_eloss,xb_new); % diimfp

    %% spectrum bulk
    data.Mat{1}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_b'));
    Rb = NSReflection(Layer(data.Mat{1}));
    Rb.theta0 = data.theta0;
    Rb.N_in = 15;
    Rb.Calculate;
    Rb.CalculateInelasticScatteringDistribution(data.theta,data.phi);
    Cb = Rb.InelasticScatteringDistribution;
%     Cb = Rb.InelasticScatteringDistribution/Rb.InelasticScatteringDistribution(1);
    convs_bulk = Rb.CalculateEnergyConvolutions;
    signal_b = sum(convs_bulk*diag(Cb),2);

    %% spectrum surface
    data.Mat{2}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_s'));
    Rs = NSReflection(Layer(data.Mat{2}));
    Rs.theta0 = data.theta0;
    Rs.N_in = 15;
    Rs.Calculate;
    convs_surf = Rs.CalculateEnergyConvolutions;
    Cs = poisspdf(0:Rs.N_in,sep);
%     Cs = Cs/Cs(1);
    signal_s = sum(convs_surf*diag(Cs),2);

    %% gauss convolution   
    signal_bs = conv_my(signal_b,signal_s,data.dE);
%     signal_sbs = conv_my(signal_s,signal_bs,data.dE);
    signal_sbs = signal_bs/max(signal_bs);
    signal_sbs_full = [signal_sbs; zeros(length(data.E0+data.dE:data.dE:data.E0+10*data.sigma_C),1)];
    signal_sbs_gauss = conv_my(signal_sbs_full,data.Gauss_C,data.dE,'same')/data.dE;
    res = signal_sbs_gauss + data.Gauss_H'*o.H + data.Gauss_D'*o.D;
    y = interp1(data.final_mesh,res*data.exp_area,xdata); 
    y(isnan(y)) = eps;
end


function [val, gradient] = fit_func_nlopt(os)
    eval_num = eval_num + 1;
    o = vecToStruct(os);

    %% x_in        
%     [x_in_b, x_in_s, int_over_depth_sigma_surf] = crosssection(o);
    dsep = dsep_Tung_2(o, data.E0, 12, false, 60);
    diimfp = ndiimfp(o, data.E0, 12, false);

    sep = trapz(o.eloss,dsep); %*h2ev;

    xs_new = interp1(o.eloss,dsep,data.mesh_eloss);
    xs_new(isnan(xs_new)) = 0;
    x_in_s = xs_new./trapz(data.mesh_eloss,xs_new); % dsep

    xb_new = interp1(o.eloss,diimfp,data.mesh_eloss);
    xb_new(isnan(xb_new)) = 0;
    x_in_b = xb_new./trapz(data.mesh_eloss,xb_new); % diimfp

    %% spectrum bulk
    data.Mat{1}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_b'));
    Rb = NSReflection(Layer(data.Mat{1}));
    Rb.theta0 = data.theta0;
    Rb.N_in = 15;
    Rb.Calculate;
    Rb.CalculateInelasticScatteringDistribution(data.theta,data.phi);
    Cb = Rb.InelasticScatteringDistribution;
%     Cb = Rb.InelasticScatteringDistribution/Rb.InelasticScatteringDistribution(1);
    convs_bulk = Rb.CalculateEnergyConvolutions;
    signal_b = sum(convs_bulk*diag(Cb),2);

    %% spectrum surface
    data.Mat{2}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_s'));
    Rs = NSReflection(Layer(data.Mat{2}));
    Rs.theta0 = data.theta0;
    Rs.N_in = 15;
    Rs.Calculate;
    convs_surf = Rs.CalculateEnergyConvolutions;
    Cs = poisspdf(0:Rs.N_in,sep);
%     Cs = Cs/Cs(1);
    signal_s = sum(convs_surf*diag(Cs),2);

    %% gauss convolution   
    signal_bs = conv_my(signal_b,signal_s,data.dE);
%     signal_sbs = conv_my(signal_s,signal_bs,data.dE);
    signal_sbs = signal_bs/max(signal_bs);
    signal_sbs_full = [signal_sbs; zeros(length(data.E0+data.dE:data.dE:data.E0+10*data.sigma_C),1)];
    signal_sbs_gauss = conv_my(signal_sbs_full,data.Gauss_C,data.dE,'same')/data.dE;
    res = signal_sbs_gauss + data.Gauss_H'*o.H + data.Gauss_D'*o.D;
    y = interp1(data.final_mesh,res*data.exp_area,data.x_exp);   
    val = sum((data.y_exp - y).^2);
    if mod(eval_num,1) == 0
        disp(['Number of function evaluations:', num2str(eval_num), ' chisq:', num2str(val)]);
    end
    if (nargout > 1)
        gradient = [0, 0.5 / val];
    end
end

function [val, gradient] = aconstraint(x)
    o = vecToStruct(x);
    o_au = convert2au(o);
    
    if strcmp(o.model,'Drude')
        [eloss,elf_henke] = mopt(o,xraypath);
        ind = bsxfun(@gt,eloss,100);
        bsum_henke = 1/(2*pi^2)*trapz(eloss(ind)/h2ev,bsxfun(@times,eloss(ind)/h2ev,elf_henke(ind)));
        cf = 4*pi*(o.ne*a0^3 - bsum_henke) / sum(o_au.A);
%         else
%             cf = o.ne * wpc / sum(o_au.A);
%         end        
    elseif strcmp(o.model,'DrudeLindhard')
        cf = (1-1/o.n_refrac^2) / sum(o_au.A);
    else
        cf = 1 / sum(o_au.A);
    end
	val = abs(cf-1);
%     if mod(eval_num,100) == 0
%         disp(['A factor:', num2str(val)]);
%     end
    if (nargout > 1)
        gradient = [0, 0.5 / val];
    end   
end

function val = fit_func_fmincon(os)

    o = vecToStruct(os);

    %% x_in        
    [x_in_b, x_in_s, int_over_depth_sigma_surf] = crosssection(o);

    %% spectrum bulk
    data.Mat{1}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_b'));
    Rb = NSReflection(Layer(data.Mat{1}));
    Rb.theta0 = data.theta0;
    Rb.N_in = 15;
    Rb.Calculate;
    Rb.CalculateInelasticScatteringDistribution(data.theta,data.phi);
    Cb = Rb.InelasticScatteringDistribution;
%     Cb = Rb.InelasticScatteringDistribution/Rb.InelasticScatteringDistribution(1);
    convs_bulk = Rb.CalculateEnergyConvolutions;
    signal_b = sum(convs_bulk*diag(Cb),2);

    %% spectrum surface
    data.Mat{2}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_s'));
    Rs = NSReflection(Layer(data.Mat{2}));
    Rs.theta0 = data.theta0;
    Rs.N_in = 15;
    Rs.Calculate;
    convs_surf = Rs.CalculateEnergyConvolutions;
    Cs = poisspdf(0:Rs.N_in,int_over_depth_sigma_surf);
%     Cs = Cs/Cs(1);
    signal_s = sum(convs_surf*diag(Cs),2);

    %% gauss convolution   
    signal_bs = conv_my(signal_b,signal_s,data.dE);
    signal_sbs = conv_my(signal_s,signal_bs,data.dE);
    signal_sbs = signal_sbs/max(signal_sbs);
    signal_sbs_full = [signal_sbs; zeros(length(data.E0+data.dE:data.dE:data.E0+10*data.sigma_C),1)];
    signal_sbs_gauss = conv_my(signal_sbs_full,data.Gauss_C,data.dE,'same')/data.dE;
    res = signal_sbs_gauss + data.Gauss_H'*data.int_H + data.Gauss_D'*data.int_D;
    y = interp1(data.final_mesh,res*data.exp_area,data.x_exp);   
    val = sum((data.y_exp - y).^2);
end

end
