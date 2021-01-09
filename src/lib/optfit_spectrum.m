function fit_result = optfit_spectrum(data,osc,xraypath)

eval_num = 0;

%% constraints
osc_min.A = ones(size(osc.A))*1e-10;
osc_min.G = ones(size(osc.G))*0.02; 
osc_min.Om = ones(size(osc.Om))*osc.Om(1);

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
osc_max.Om = ones(size(osc.Om))*1000;

if osc.H ~= 0
    osc_min.H = eps;
    osc_max.H = 0.05;
end
if osc.D ~= 0
    osc_min.D = eps;
    osc_max.D = 0.05;
end

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
opt.maxeval = 10000;
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
%{
figure;
plot(data.E0 - data.x_exp,data.y_exp,'DisplayName','Experiment','Marker','o','LineWidth',1)
hold on
plot(data.E0 - data.x_exp,fit_func(x_res),'DisplayName','Summary signal','LineWidth',2)
legend

xlabel('Kinetic energy, eV')
ylabel('Intensity, rel.un.')
set(findall(gcf,'-property','FontSize'),'FontSize',25)
%}

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
    if osc.H ~= 0 && osc.D == 0
        v = [s.A, s.G, s.Om, s.H];
    elseif osc.D ~= 0 && osc.H == 0
        v = [s.A, s.G, s.Om, s.D];
    else
        v = [s.A, s.G, s.Om, s.H, s.D];
    end
end

function s = vecToStruct(v)
    s = osc;
    nA = length(osc.A);

    s.A = v(1:nA);
    s.G = v((1:nA)+nA);
    s.Om = v((1:nA)+nA*2);
    if osc.H ~= 0 && osc.D == 0
        s.H = v(end);
    elseif osc.D ~= 0 && osc.H == 0
        s.D = v(end);
    else
        s.H = v(end-1);
        s.D = v(end);
    end
end

function [x_in_b, x_in_s, sep] = crosssection(o)
    depth = -5:5;
    clear_bulk = zeros(length(o.eloss), length(depth));
    reduced_bulk = zeros(length(o.eloss), length(depth));
    surface = zeros(length(o.eloss), length(depth));
    E0 = data.E0;
    theta = data.theta;

    parfor r = 1:length(depth)
        disp(r)
        [clear_bulk(:,r),surface(:,r),reduced_bulk(:,r)] = ndiimfp_Li(o,E0,depth(r),theta,9);
    end
    ind = depth >= 0;
    bulk = clear_bulk(:,ind) + reduced_bulk(:,ind);
    bulk = trapz(depth(ind),bulk,2);
    surface = trapz(depth,surface,2);

    sep = trapz(o.eloss,surface); 

    xs_new = interp1(o.eloss,surface,data.mesh_eloss);
    xs_new(isnan(xs_new)) = 0;
    x_in_s = xs_new./trapz(data.mesh_eloss,xs_new); % dsep

    xb_new = interp1(o.eloss,bulk,data.mesh_eloss);
    xb_new(isnan(xb_new)) = 0;
    x_in_b = xb_new./trapz(data.mesh_eloss,xb_new); % diimfp
end

function y = fit_func(os)

    o = vecToStruct(os);

    %% x_in        
    [x_in_b, x_in_s, sep] = crosssection(o);

    %% spectrum bulk
    data.Mat{1}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_b'));
    Rb = NSReflection(Layer(data.Mat{1}));
    Rb.theta0 = data.theta0;
    Rb.N_in = 10;
    Rb.Calculate;
    Rb.CalculateInelasticScatteringDistribution(data.theta,data.phi);
    Cb = Rb.InelasticScatteringDistribution/Rb.InelasticScatteringDistribution(1); % !!! very important to get the inelastic background right
    convs_bulk = Rb.CalculateEnergyConvolutions;
    signal_b = sum(convs_bulk*diag(Cb),2);

    %% spectrum surface
    data.Mat{2}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_s'));
    Rs = NSReflection(Layer(data.Mat{2}));
    Rs.theta0 = data.theta0;
    Rs.N_in = 10;
    Rs.Calculate;
    convs_surf = Rs.CalculateEnergyConvolutions;
    Cs = poisspdf(0:Rs.N_in,sep);
    signal_s = sum(convs_surf*diag(Cs),2);

    %% gauss convolution   
    signal_bs = conv_my(signal_b,signal_s,data.dE); % bulk+surface convolution
    signal_bs_full = [signal_bs; zeros(length(data.E0+data.dE:data.dE:data.E0+10*data.sigma_C),1)];
    signal_bs_full_gauss = conv_my(signal_bs_full,data.Gauss_C,data.dE);
    signal_bs_full_gauss_h = signal_bs_full_gauss + data.Gauss_H'*osc.H; %+ data.Gauss_D'*osc.D
    y = interp1(data.final_mesh,signal_bs_full_gauss_h,data.x_exp);  
    y(isnan(y)) = eps;
end


function [val, gradient] = fit_func_nlopt(os)
    
    eval_num = eval_num + 1;
    o = vecToStruct(os);

    %% x_in        
    [x_in_b, x_in_s, sep] = crosssection(o);

    %% spectrum bulk
    data.Mat{1}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_b'));
    Rb = NSReflection(Layer(data.Mat{1}));
    Rb.theta0 = data.theta0;
    Rb.N_in = 10;
    Rb.Calculate;
    Rb.CalculateInelasticScatteringDistribution(data.theta,data.phi);
    Cb = Rb.InelasticScatteringDistribution/Rb.InelasticScatteringDistribution(1); % !!! very important to get the inelastic background right
    convs_bulk = Rb.CalculateEnergyConvolutions;
    signal_b = sum(convs_bulk*diag(Cb),2);

    %% spectrum surface
    data.Mat{2}.SetManualDIIMFP(flipud((data.E0-data.mesh_eloss)'),flipud(x_in_s'));
    Rs = NSReflection(Layer(data.Mat{2}));
    Rs.theta0 = data.theta0;
    Rs.N_in = 10;
    Rs.Calculate;
    convs_surf = Rs.CalculateEnergyConvolutions;
    Cs = poisspdf(0:Rs.N_in,sep);
    signal_s = sum(convs_surf*diag(Cs),2);

    %% gauss convolution   
    signal_bs = conv_my(signal_b,signal_s,data.dE); % bulk+surface convolution
    signal_bs_full = [signal_bs; zeros(length(data.E0+data.dE:data.dE:data.E0+10*data.sigma_C),1)];
    signal_bs_full_gauss = conv_my(signal_bs_full,data.Gauss_C,data.dE);
    signal_bs_full_gauss_h = signal_bs_full_gauss + data.Gauss_H'*osc.H; %+ data.Gauss_D'*osc.D
    y = interp1(data.final_mesh,signal_bs_full_gauss_h,data.x_exp); 
    y(isnan(y)) = eps;
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

end
