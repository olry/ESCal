function fit_result = optfit(x_exp,y_exp,osc,fitgoal,E0,varargin)

if nargin<5 && strcmp(fitgoal,'ndiimfp')
    error('Set incident energy E0')
end

osc_min.A = zeros(size(osc.A));
osc_min.G = ones(size(osc.G))*0.02; 
osc_min.Om = zeros(size(osc.A));

osc_max.A = ones(size(osc.A))*Inf;
osc_max.G = ones(size(osc.A))*100;
osc_max.Om = ones(size(osc.A))*x_exp(end);

lb = structToVec(osc_min);
ub = structToVec(osc_max);

%% LSQ fitting (local minimum search)
%{
options = optimoptions('lsqcurvefit','OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16,'PlotFcn',@optimplotx,'Display','iter-detailed','UseParallel',true);
options.FunctionTolerance = 1.000000e-12;
options.OptimalityTolerance = 1.000000e-12;
options.MaxFunctionEvaluations = 7000;
options.StepTolerance = 1e-12;
x_res = lsqcurvefit(@fit_func, structToVec(osc), x_exp, y_exp, lb, ub, options);
%}
%% Global

rng default % For reproducibility
% %gs = GlobalSearch;
% 
problem = createOptimProblem('fmincon',...
   'objective',@fit_func_2,...
   'x0',structToVec(osc),'options',...
   optimoptions('fmincon','Algorithm','sqp','Display','iter','UseParallel',true));
problem.lb = lb;
problem.ub = ub;
problem.options.MaxIterations = 1e5;
problem.options.OptimalityTolerance = 1e-15;
problem.options.StepTolerance = 1e-15;
problem.options.MaxFunctionEvaluations = 1e5;
% % [x_res] = run(gs,problem);
% 
x_res = fmincon(problem);

an = vecToStruct(x_res);
disp(an);

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

an = scaling(an);
fit_result = an;
an.eloss = [eps:0.1:10000]';
an_au = convert2au(an);
bsum = 1/(2*pi^2)*trapz(an_au.eloss,bsxfun(@times,an_au.eloss,eps_sum(an)));
psum = 2/pi*trapz(an.eloss,bsxfun(@rdivide,eps_sum(an),an.eloss)) + 1/an.n_refrac^2;

disp(['P-sum rule: ',num2str(psum)]);
disp(['Bethe sum rule: ',num2str(bsum), ' electron density = ',num2str(osc.ne*a0^3), '(a.u.^-3)']);
disp(['Sum of A: ',num2str(sum(an_au.A)),'(1-1/n2) = ', num2str(abs(1-(1/osc.n_refrac)^2)), ' 4piN = ',num2str(4*pi*osc.ne*a0^3), '(a.u.^-3)']);

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

function osc_scaled = scaling(o)
    o = convert2au(o);
    switch o.model
        case 'Drude'
            o.A = o.A/sum(o.A)*4*pi*(o.ne*a0^3);
            eps1 = eps_real_Drude(o);
            w_ScalingFactor = sqrt((eps1(1) - o.beps)/(o.n_refrac^2 - o.beps));
            o.Om = o.Om / w_ScalingFactor;
        case 'DrudeLindhard'
            o.A = o.A/sum(o.A)*(1-1/o.n_refrac^2);
            BetheSum = sum((pi/2)*o.A.*o.Om.^2);
            BetheValue = 2*pi*pi*o.ne*a0^3;
            w_ScalingFactor = sqrt(BetheSum / BetheValue);
            o.Om = o.Om / w_ScalingFactor;
        case 'Mermin'
            o.A = o.A/sum(o.A); %*(1-1/o.n_refrac^2);
            BetheSum = sum((pi/2)*o.A.*o.Om.^2);
            BetheValue = 2*pi*pi*o.ne*a0^3;
            w_ScalingFactor = sqrt(BetheSum / BetheValue);
            o.Om = o.Om / w_ScalingFactor;
    end
    osc_scaled = convert2ru(o);
end

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

function y = fit_func(x,xdata)
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
            y = interp1(o.eloss,diimfp,xdata);
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
