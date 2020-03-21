function fit_result = optfit_ndiimfp(x_exp,y_exp,E0,osc)

osc_min.A = zeros(size(osc.A));
osc_min.G = ones(size(osc.G))*0.02; 
osc_min.Om = zeros(size(osc.A));

osc_max.A = Inf(size(osc.A));
osc_max.G = Inf(size(osc.A));
osc_max.Om = Inf(size(osc.A));

lb = structToVec(osc_min);
ub = structToVec(osc_max);

%% LSQ fitting (local minimum search)

% options = optimoptions('lsqcurvefit','OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16,'PlotFcn',@optimplotx,'Display','iter-detailed','UseParallel',true);
% options.FunctionTolerance = 1.000000e-12;
% options.OptimalityTolerance = 1.000000e-12;
% options.MaxFunctionEvaluations = 7000;
% options.StepTolerance = 1e-12;
% pars = structToVec(osc);
% x_res = lsqcurvefit(@fit_func, pars, x_exp, y_exp, lb, ub, options);

%% Global

rng default 
problem = createOptimProblem('fmincon',...
   'objective',@fit_func_2,...
   'x0',structToVec(osc),'options',...
   optimoptions('fmincon','Algorithm','sqp','Display','iter','UseParallel',true));
problem.lb = lb;
problem.ub = ub;
problem.options.MaxIterations = 1e4;
problem.options.OptimalityTolerance = 1e-12;
problem.options.StepTolerance = 1e-12;
problem.options.MaxFunctionEvaluations = 1e6;

x_res = fmincon(problem);

an = vecToStruct(x_res);
disp(an);

%% Plotting

h = figure;
plot(x_exp,y_exp,'DisplayName','Experimental nDIIMFP','Marker','o','LineWidth',1)
hold on
plot(x_exp,fit_func(x_res,x_exp),'DisplayName','This work','LineWidth',2);

ylabel('nDIIMFP');
xlabel('Energy loss $\omega$, eV');

title(osc.name);
set(findall(gcf,'-property','FontSize'),'FontSize',24)
xlim([0 150]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
txt = [osc.name, '_fit'];
print(h,txt,'-dpdf','-r0')
txt = [osc.name,'_fit.fig'];
savefig(txt)

%% Sum-rules

an = scaling_drude(an);
fit_result = an;

an_au = convert2au(an);
bsum = 1/(2*pi^2)*trapz(an_au.eloss,bsxfun(@times,an_au.eloss,eps_sum(an)));
psum = 2/pi*trapz(an_au.eloss,bsxfun(@rdivide,eps_sum(an),an_au.eloss));

disp(['P-sum rule: ',num2str(psum)]);
disp(['Bethe sum rule: ',num2str(bsum), ' electron density = ',num2str(osc.ne*a0^3), '(a.u.^-3)']);
disp(['Sum of A: ',num2str(sum(an_au.A)), ' 4piN = ',num2str(4*pi*osc.ne*a0^3), '(a.u.^-3)']);

fsum = 1/(2*pi^2*(osc.na*a0^3))*trapz(an_au.eloss,bsxfun(@times,an_au.eloss,eps_sum(an)));
disp(['F-sum rule:',num2str(fsum), ' EAN = ', num2str(osc.Z)]);



%% Functions

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

function osc_scaled = scaling_dl(o)
    osc_scaled = o;
    osc_scaled.A = o.A/sum(o.A)*(1-(1/osc.n_refrac)^2);
    elf0 = eps_sum(o);
    elf0(isnan(elf0)) = 0;
    o_au = convert2au(o);
    N = 1/(2*pi^2)*trapz(o_au.eloss',o_au.eloss.*elf0(:,1)');        
    BetheSum = sum((pi/2)*o_au.A.*o_au.Om.^2);
    BetheValue = 2*pi*pi*ne;
    w_ScalingFactor = sqrt(BetheSum / BetheValue);
    osc_scaled.Om = o.Om / w_ScalingFactor;
end

function osc_scaled = scaling_drude(o)
    o = convert2au(o);
    o.A = o.A/sum(o.A)*4*pi*(o.ne*a0^3);
    eps1 = eps_real_Drude(o);
    w_ScalingFactor = sqrt((eps1(1) - o.beps)/(o.n_refrac^2 - o.beps));
    o.Om = o.Om / w_ScalingFactor;
    osc_scaled = convert2ru(o);
end

function eps = eps_real_Drude(o)
    %o = convert2au(o);
    w = o.eloss;
    q = o.qtran;
    eps_re = o.beps;
    for j=1:length(o.A)
        [epsDrud_re, epsDrud_im] = Drude(q,w,o.Om(j),o.G(j),o.alpha,o.Ef);
        eps_re = eps_re - o.A(j)*epsDrud_re;
    end
    eps = eps_re;
end

function y = fit_func(x,xdata)
    o = vecToStruct(x);
    o.eloss = [0:.5:200 250:10:E0];
    o = scaling_drude(o);
    diimfp = ndiimfp(o,E0,10);
    diimfp(1) = eps;
    y = interp1(o.eloss,diimfp,xdata);
end

function err = fit_func_2(x)
    o = vecToStruct(x);
    o.eloss = [0:.5:200 250:10:E0];
    o = scaling_drude(o);
    diimfp = ndiimfp(o,E0,10);
    diimfp(1) = eps;
    y = interp1(o.eloss,diimfp,x_exp);
    weight = 10;
    err = trapz(x_exp, weight*abs(y_exp-y));
end

end
