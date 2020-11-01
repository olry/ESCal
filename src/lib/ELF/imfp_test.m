clear;

E = [20:100 200:100:1000 2000 3000 4000];
lambda = zeros(size(E));

%{
osc.model = 'Drude';
osc.A = [245 600];
osc.G = [1 1];
osc.Om = [0 76];
osc.alpha = 1;
osc.Ef = 10.9;
osc.egap = 7;
osc.vb = 10;
osc.beps = 1;
osc.qtran = 0;
%}
%{
osc.model = 'DrudeLindhard';
osc.A = [0.9 0.1];
osc.G = [1 1];
osc.Om = [15 80];
osc.alpha = 1;
osc.Ef = 10.9;
osc.egap = 7;
osc.vb = 10;
osc.beps = 1;
osc.qtran = 0;
%}
% {
osc.model = 'Mermin';
osc.A = [0.9 0.1];
osc.G = [1 1];
osc.Om = [15 80];
osc.alpha = 1;
osc.Ef = 10.9;
osc.egap = 7;
osc.vb = 10;
osc.beps = 1;
osc.qtran = 0.01;
%}

for i = 1:length(E)
%     lambda(i) = imfp(osc,E(i));
    lambda(i) = imfp(osc,E(i),'metal',true);
end

% figure
semilogx(E,lambda,'LineWidth',2,'DisplayName',osc.model)
% xlim([20 4000])
% ylim([0 60])


