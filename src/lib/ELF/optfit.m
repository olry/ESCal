function fit_result = optfit(x_exp,y_exp,osc,fitgoal,E0,xraypath,varargin)

eval_num = 0;

osc_min.A = ones(size(osc.A))*1e-10;
osc_min.G = ones(size(osc.G))*0.25/h2ev; 
osc_min.Om = ones(size(osc.Om))*osc.egap;

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

%% LSQ fitting (local minimum search)
%{
rng default;

options = optimoptions('lsqcurvefit','UseParallel',true,'Display','iter');
% options.FunctionTolerance = 1.000000e-12;
% options.OptimalityTolerance = 1.000000e-12;
% options.MaxFunctionEvaluations = 1e4;
% options.StepTolerance = 1e-12;
% options.MaxIterations = 1e4;
x_res = lsqcurvefit(@fit_func_elf, structToVec(osc), x_exp, y_exp, lb, ub, options);
an = vecToStruct(x_res);
an = scaling(an);
% an = scaling_ohne_henke(an);
%}
%% NLopt
% {
opt.algorithm = NLOPT_LN_COBYLA;
opt.lower_bounds = lb;
opt.upper_bounds = ub;
opt.maxeval = 1e6;
opt.min_objective = @fit_func_nlopt;
if osc.henke
    if strcmp(osc.model,'Drude')
       opt.fc = { (@(x) aconstraint_henke(x)) };
    else
        opt.fc = { (@(x) wconstraint_henke(x)); };
    end
else
    if strcmp(osc.model,'Drude')
       opt.fc = { (@(x) aconstraint(x)) };
    else
        opt.fc = { (@(x) wconstraint(x)) };
    end
end
opt.fc_tol = 1e-8; 
opt.xtol_rel = 1e-10;
[x_res] = nlopt_optimize(opt, structToVec(osc));
an = vecToStruct(x_res);
%}
%% Fmincon solver
%{
rng default % For reproducibility

if globalSearch
%     ms = MultiStart('FunctionTolerance',2e-4,'UseParallel',true);
    ms = MultiStart('UseParallel',true);
    gs = GlobalSearch(ms);
    problem = createOptimProblem('fmincon',...
    'objective',@fit_func_2,...
    'x0',structToVec(osc),'options',...
    optimoptions('fmincon','Algorithm','sqp','Display','iter','UseParallel',true));
    problem.lb = lb;
    problem.ub = ub;
    x_res = run(gs,problem);
%     x_res = run(ms,problem,20);
else
    problem = createOptimProblem('fmincon',...
    'objective',@fit_func_2,...
    'x0',structToVec(osc),'options',...
    optimoptions('fmincon','Algorithm','sqp','Display','iter','UseParallel',true));
    problem.lb = lb;
    problem.ub = ub;
%     problem.options.MaxIterations = 1e6;
%     problem.options.OptimalityTolerance = 1e-15;
%     problem.options.StepTolerance = 1e-15;
%     problem.options.MaxFunctionEvaluations = 1e6;
    x_res = fmincon(problem);
end
%}

disp(an);
fit_result = an;

%% Plotting
h = figure;
plot(x_exp,y_exp,'DisplayName','Experimental nDIIMFP','Marker','o','LineWidth',1)
hold on
plot(x_exp,fit_func(x_res,x_exp),'DisplayName','This work','LineWidth',2);

ylabel('nDIIMFP');
xlabel('Energy loss $\omega$ (eV)');
txt = [osc.name,' ',osc.model];
title(txt);
set(findall(gcf,'-property','FontSize'),'FontSize',24)
xlim([0 100]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
txt = [osc.name, '_fit_',osc.model];
print(h,txt,'-dpdf','-r0')
txt = [osc.name,'_fit_',osc.model,'.fig'];
savefig(txt)

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

%% Functions

function v = structToVec(s)
    if strcmp( osc.model,'MerminLL')
        v = [s.A, s.G, s.Om, s.u];
    else
        v = [s.A, s.G, s.Om];
    end
end

function s = vecToStruct(v)   
    s = osc;    
    nA = length(osc.A);    
    s.A = v(1:nA);
    s.G = v((1:nA)+nA);
    s.Om = v((1:nA)+nA*2);
    if strcmp( osc.model,'MerminLL')
        s.u = v(end);
    end
end

function y = fit_func(x,xdata)
    o = vecToStruct(x);
%     o = scaling(o);
    switch fitgoal
        case 'elf'           
            [eloss,elf] = mopt(o,xraypath,true);
            y = interp1(eloss,elf,xdata);
        case 'ndiimfp'
            diimfp = ndiimfp(o,E0,10,true,false,xraypath);
            diimfp(1) = eps;
            y = interp1(o.eloss,diimfp,xdata);
    end
end

function [val, gradient] = fit_func_nlopt(x)
    eval_num = eval_num + 1;
    o = vecToStruct(x);
    switch fitgoal
        case 'elf'
            [eloss,elf] = mopt(o,xraypath,true);
            res = interp1(eloss,elf,x_exp);
        case 'ndiimfp'
            diimfp = ndiimfp(o,E0,10,true,false,xraypath);
            diimfp(1) = eps;
            res = interp1(o.eloss,diimfp,x_exp);
    end
    val = sum((y_exp - res).^2);
    if mod(eval_num,100) == 0
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
        cf = o.ne * wpc / sum(o_au.A);
    elseif strcmp(o.model,'DrudeLindhard')
        cf = (1-1/o.n_refrac^2) / sum(o_au.A);
    else
        cf = 1 / sum(o_au.A);
    end
	val = abs(cf-1);
    if mod(eval_num,100) == 0
        disp(['A factor:', num2str(val)]);
    end
    if (nargout > 1)
        gradient = [0, 0.5 / val];
    end   
end

function [val, gradient] = wconstraint(x)
    o = vecToStruct(x);
    o_au = convert2au(o);
    if strcmp(o.model,'Drude')
        w_ScalingFactor = 1;
    else
        BetheSum = sum((pi/2)*o_au.A.*o_au.Om.^2);
        BetheValue = 2*pi*pi*(o.ne*a0^3); 
        w_ScalingFactor = sqrt(BetheSum / BetheValue);
    end
	val = abs(w_ScalingFactor-1);
    if mod(eval_num,100) == 0
        disp(['W factor:', num2str(val)]);
    end
    if (nargout > 1)
        gradient = [0, 0.5 / val];
    end
end

function [val, gradient] = aconstraint_henke(x)
    o = vecToStruct(x);
    o_au = convert2au(o);
    
    if strcmp(o.model,'Drude')
        [eloss,elf_henke] = mopt(o,xraypath);
        ind = bsxfun(@gt,eloss,100);
        bsum_henke = 1/(2*pi^2)*trapz(eloss(ind)/h2ev,bsxfun(@times,eloss(ind)/h2ev,elf_henke(ind)));
        cf = (o.ne*a0^3 - bsum_henke) * 4*pi / sum(o_au.A);
    elseif strcmp(o.model,'DrudeLindhard')
        cf = (1-1/o.n_refrac^2) / sum(o_au.A);
    else
        cf = 1 / sum(o_au.A);
    end
	val = abs(cf-1);
    if mod(eval_num,100) == 0
        disp(['A factor:', num2str(val)]);
    end
    if (nargout > 1)
        gradient = [0, 0.5 / val];
    end   
end

function [val, gradient] = wconstraint_henke(x)
    o = vecToStruct(x);
    o_au = convert2au(o);
    if strcmp(o.model,'Drude')
        w_ScalingFactor = 1;
    else
        switch fitgoal
            case 'elf'
                ind = bsxfun(@le,x_exp,100);
                bsum_henke = 1/(2*pi^2)*trapz(x_exp(ind)/h2ev,bsxfun(@times,x_exp(ind)/h2ev,y_exp(ind)));
                BetheSum = sum((pi/2)*o_au.A.*o_au.Om.^2);
                BetheValue = 2*pi*pi*bsum_henke; 
            case 'ndiimfp'
                [eloss,elf_henke] = mopt(o,xraypath);
                ind = bsxfun(@gt,eloss,100);
                bsum_henke = 1/(2*pi^2)*trapz(eloss(ind)/h2ev,bsxfun(@times,eloss(ind)/h2ev,elf_henke(ind)));
                BetheSum = sum((pi/2)*o_au.A.*o_au.Om.^2);
                BetheValue = 2*pi*pi*(o.ne*a0^3 - bsum_henke); 
        end      
        w_ScalingFactor = sqrt(BetheSum / BetheValue);
    end
	val = abs(w_ScalingFactor-1);
    if mod(eval_num,100) == 0
        disp(['W factor:', num2str(val)]);
    end
    if (nargout > 1)
        gradient = [0, 0.5 / val];
    end
end

function err = fit_func_2(x)
    o = vecToStruct(x);
    o = scaling(o);
    switch fitgoal
        case 'elf'           
            elf = eps_sum(o);
            elf(1) = eps;
            y = elf;
        case 'ndiimfp'
            diimfp = ndiimfp(o,E0,10);
            diimfp(1) = eps;
            y = interp1(o.eloss,diimfp,x_exp);
    end
    weight = 1;
    err = trapz(x_exp, weight*abs(y_exp - y));
end

end
