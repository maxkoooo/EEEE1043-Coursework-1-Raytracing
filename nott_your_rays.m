classdef nott_your_rays
    properties
        f {mustBeNumeric}           % signal frequency
        tx {mustBeNumeric}          % transmitter coordinates
        rx {mustBeNumeric}          % receiver coordinates
        rxf {mustBeNumeric}         % final receiver coordinates
        twy {mustBeNumeric}         % y-level of top wall
        bwy {mustBeNumeric}         % y-level of bottom wall
        resolution {mustBeNumeric}  % resolution of receiver movement
        epsilon1 {mustBeNumeric}    % permittivity of air
        epsilon2 {mustBeNumeric}    % permittivity of reflection surface
        n {mustBeNumeric}           % number of orders of reflections to calculate
        plotMag {mustBeNumericOrLogical}     % switch to plot magnitude
    end

    methods
        function [Eddb, Esdb, rxv] = raytrace(obj)
            %% preamble
            c = 3*10^8;     % speed of light
            beta = 2*pi*obj.f/c;    % wave number
            rxv = obj.rx(1):obj.resolution:obj.rxf;    % move receiver from p to q (12 meters)
            rxn = ((obj.rxf-obj.rx(1))/obj.resolution) + 1;
            relperm = obj.epsilon2/obj.epsilon1;   % relative permittivity (dielectric constant) between the wall and the air
            index = obj.n+1;    % for indexing the transmitter image array, where index of 1 is set to the original coordinates of the transmitter
            anteff = (c/(4*pi*obj.f));  % antenna effect (shorter antenna needed to capture higher frequency, total intensity of received signal is lesser)

            %% array pre-allocation
            % pre-allocated for speed
            % columns vary with rxv, rows vary with order of reflection
            % page 1 for top-wall-first reflections, page 2 for bottom-wall-first reflections

            ity = zeros(1, index, 2);   % images of ty
            diy = zeros(1, obj.n, 2);
            di = zeros(obj.n, rxn, 2);
            refangle = zeros(obj.n, rxn, 2);
            transtheta = zeros(obj.n, rxn, 2);
            gamma = zeros(obj.n, rxn, 2);
            Ei = zeros(obj.n, rxn, 2);
            Eidb = zeros(obj.n, rxn, 2);
            Eis = zeros(obj.n, rxn);
            Es = zeros(obj.n, rxn);
            Esdb = zeros(obj.n, rxn);
            
            ity(1, 1, :) = obj.tx(2);  % initialise first element to ty
            
            %% image calculation
            % calculate y axis coordinates for all images
            for i = 2:index % start from 2 as 1 is indexed as the original coordinates of ty
                m = mod(i, 2);   % boolean true for odd number of i
                h = not(m);      % boolean true for even number of i
                
                % first reflection on the top wall
                ity(1, i, 1) = ity(1, (i-1), 1) + (2*(obj.twy-ity(1, (i-1), 1))*h) + (2*(obj.bwy-ity(1, (i-1), 1))*m);
                
                % first reflection on the bottom wall
                ity(1, i, 2) = ity(1, (i-1), 2) + (2*(obj.twy-ity(1, (i-1), 2))*m) + (2*(obj.bwy-ity(1, (i-1), 2))*h);
            end
            
            %% direct ray calculation
            ddirect = sqrt(((rxv-obj.tx(1)).^2) + ((obj.rx(2)-obj.tx(2))^2) + ((obj.rx(3)-obj.tx(3))^2));   % distance between transmitter and receiver
            Ed = (1./ddirect).*exp(-1j*beta*ddirect).*anteff;    % electric field strength of direct ray
            Eddb = 20*log10(abs(Ed));   % electric field strength of direct ray in dB
            
            %% reflected rays calculation
            % each column of the arrays are made to vary with rxv
            % each row of the array varies with the order of reflection
            for i = 1:obj.n % iterates through the order of reflection
                for k = 1:2 % iterates through the pages of the array, 1 for top-wall-first reflections, and 2 for bottom-wall-first reflections
                    diy(1, i, k) = ity(1, (i+1), k)-obj.rx(2);   % y-axis difference between top-first reflection transmitter image and receiver
                
                    % distance between image and receiver
                    di(i, :, k) = sqrt((rxv.^2) + (diy(1, i, k)^2) + ((obj.rx(3)-obj.tx(3))^2)); % distance between image of transmitter and receiver through n number of orders of reflections, with first reflection along the top wall
                
                    % angle of reflection calculated using image distance from receiver (hypothenuse) and x-axis distance of the image to the receiver
                    % property of similar triangles means angles are the same
                    % angle of incidence is equal to the angle of reflection
                    % angles for a particular order of reflection at a particular distance are same for each reflection
                    refangle(i, :, k) = asin(rxv./di(i, :, k));
                
                    % angle of transmission, found using snell's law
                    transtheta(i, :, k) = asin(sin(refangle(i, :, k))/sqrt(relperm));    
                
                    % reflection coefficient per reflection
                    % gamma for each reflection in the same order and distance is equivalent
                    gamma(i, :, k) = (cos(refangle(i, :, k)) - (sqrt(relperm).*cos(transtheta(i, :, k))))./(cos(refangle(i, :, k)) + sqrt(relperm).*cos(transtheta(i, :, k)));
                
                    % electric field strength per order of reflection
                    % gamma is raised to the power of the order of reflection, i.e. how many times the ray reflects off the wall
                    Ei(i, :, k) = (gamma(i, :, k).^i).*(1./di(i, :, k)).*exp(-1j*beta.*di(i, :, k)).*anteff;   
                
                    % electric field strength per order of reflection plotted in the dB scale
                    Eidb(i, :, k) = 20*log10(abs(Ei(i, :, k)));    
                end  
            
                % cumulative sum of electric field strength for an increasing number of orders of reflection
                if obj.plotMag == 1 % to sum up only the magnitude of the electric fields
                    if i == 1   
                        Eis(i, :) = abs(Ei(i, :, 1)) + abs(Ei(i, :, 2)); % for initialising the first element
                    else
                        Eis(i, :) = Eis((i-1), :) + abs(Ei(i, :, 1)) + abs(Ei(i, :, 2));
                    end

                    % total sum of the electric field strength of the reflections and the direct ray as seen at the receiver
                    Es(i, :) = abs(Ed) + Eis(i, :);
                else    
                    if i == 1   
                        Eis(i, :) = Ei(i, :, 1) + Ei(i, :, 2); % for initialising the first element
                    else
                        Eis(i, :) = Eis((i-1), :) + Ei(i, :, 1) + Ei(i, :, 2);
                    end

                    % total sum of the electric field strength of the reflections and the direct ray as seen at the receiver
                    Es(i, :) = Ed + Eis(i, :);
                end
            
                % total electric field strength plotted in the dB scale
                Esdb(i, :) = 20*log10(abs(Es(i, :)));
            end
        end

        function lgd = dynamicLegends(~, i)
            if i == 1
                lgd = strcat(num2str(i), 'st Order');
            elseif i == 2
                lgd = strcat(num2str(i), 'nd Order');
            elseif i == 3
                lgd = strcat(num2str(i), 'rd Order');
            else
                lgd = strcat(num2str(i), 'th Order');
            end
        end
    end
end
