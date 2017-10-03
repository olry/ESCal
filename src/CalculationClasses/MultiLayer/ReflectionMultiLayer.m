classdef ReflectionMultiLayer<BaseMultiLayer
    properties
        %          z; %�������� ��� tau �� z, ������ ������ � ���������� �����
        theta                  = 0;
        phi                    = [];
        %         delta_max              = [];
        Fm;

    end
    
    methods
        function obj = ReflectionMultiLayer (Layers, CalculationMethod)
            
            if nargin==1 || isempty(CalculationMethod);
                CalculationMethod = {'NS'};             
            end
            
            obj = obj@BaseMultiLayer(Layers, CalculationMethod);
            
            obj.dE = abs(obj.Layers(1).Material.DIIMFP_E(2) - obj.Layers(1).Material.DIIMFP_E(1));
            obj.energy_mesh_full=obj.Layers(1).Material.DIIMFP_E;
            % �������� �� ������������ E0  ��� ���� ����
            %             for i_layer = 1:N_Layer
            %
            %             end
            %             for i = 1:obj.N_layer
            %                     obj.z(i) = Layers(i).thickness;
            %                     obj.Mat(i)=Layers(i).MaterialName;
            %             end
        end
        
        function Calculate(obj)

              
            for i_layer = 1:obj.N_Layer
                switch obj.CalculationMethod{i_layer}
                    case 'NS'
                        obj.ObjectsOfLayers{i_layer} = NSTransmition(obj.Layers(i_layer));
                    case 'SLA'
                        obj.ObjectsOfLayers{i_layer} = SLATransmition(obj.Layers(i_layer));
                    case 'SA'
                        obj.ObjectsOfLayers{i_layer} = SATransmition(obj.Layers(i_layer));
                    otherwise
                        error('������� ���������� �������� ������: NS, SLA ��� SA.');
                end
            end
            obj.CalculateEachLayer;
            obj.mu_mesh = obj.ObjectsOfLayers{1}.mu_mesh; % !!!!!

            obj.Setup_energy_mesh_full;
            obj.RecalculateToEnergyDistribution;
            R = obj.FullEnergyDistribution{obj.N_Layer}.R(:,:,:,:);
%             w = sparse(1:obj.ObjectsOfLayers{1}.N, 1:obj.ObjectsOfLayers{1}.N, obj.ObjectsOfLayers{1}.mu_mesh_weights./obj.mu_mesh);
            M = size(R,4)-1;
            N_E =size(R,3);
            wT = zeros(size(obj.FullEnergyDistribution{1}.T));
            Tw = zeros(size(obj.FullEnergyDistribution{1}.T));
            for i_layer = (obj.N_Layer-1):-1:1
                Rlocal = obj.FullEnergyDistribution{i_layer}.R(:,:,:,:);
                if strcmp(obj.CalculationMethod{i_layer},'SLA')
                    % � SLA-������� ��� ���������� �������, ��������� ����
                    % �� ������� ������� ���������, � ���� �������
                    % ������������ �� ���� ������-�������, ������� �� �����
                    % ������������� � �� ����� ��������� ������� ������,
                    % ������� ��������� �������, ��������� ���� �� �������
                    % ������� ���������
                    wT = obj.FullEnergyDistribution{i_layer}.T;
                    Tw = obj.FullEnergyDistribution{i_layer}.T;
                else
                    w = sparse(1:obj.ObjectsOfLayers{i_layer}.N, 1:obj.ObjectsOfLayers{i_layer}.N, obj.ObjectsOfLayers{i_layer}.mu_mesh_weights./obj.mu_mesh);
                    for m=0:M
                        for i=1:N_E
                            wT(:,:,i,m+1) = w*obj.FullEnergyDistribution{i_layer}.T(:,:,i,m+1) + obj.FullEnergyDistribution{i_layer}.L(:,:,i,m+1);
                            Tw(:,:,i,m+1) = obj.FullEnergyDistribution{i_layer}.T(:,:,i,m+1)*w + obj.FullEnergyDistribution{i_layer}.L(:,:,i,m+1);
                        end
                    end
                end
                
                [~, K_E]=min(abs(obj.energy_mesh_full-obj.ObjectsOfLayers{1}.Material.E0));
                K_E = length(obj.energy_mesh_full)-K_E;%find(obj.energy_mesh_full==obj.ObjectsOfLayers{1}.Material.E0); % ����� ��������� ����� E0
 
%                 K_E = length(obj.energy_mesh_full)-find(obj.energy_mesh_full==obj.ObjectsOfLayers{1}.Material.E0); % ����� ��������� ����� E0
                R = obj.DoubleLayer(Rlocal,Tw,R,wT,K_E);
            end
            obj.Fm = R;
            obj.IsCalculated = true;
        end
        
        function CalculateEnergyDistribution(obj,theta,phi)
            if nargin < 3; phi = 0; end; %??????
            obj.CalculateEnergyDistribution_E0(theta,phi);
            convGauss (obj,obj.sigma_gauss);
            convLDS (obj,obj.sigma_LDS, obj.alpha_LDS);
        end
        
    end % methods
end % classdef