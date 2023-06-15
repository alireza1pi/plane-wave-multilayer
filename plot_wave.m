function EE = plot_wave(thickness, d , kd, theta, mag_plus, mag_minus, minY, maxY, minC, maxC, Nx, Ny, nnLayer, omega,selectedDirection, hObject, handles)

thickness = thickness(1:nnLayer);
total_Length = sum(thickness);
total_Width = maxY - minY;

kd = kd(1:nnLayer);
theta = theta(1:nnLayer);
kz = kd.*cosd(theta);
ky = kd.*sind(theta);



for ii = 1: length(thickness)
    
    if ii ==1
        x = (0 : total_Length/Nx :thickness(ii));
    else
        x = (d(ii-1) : total_Length/Nx : d(ii-1)+thickness(ii));
    end
    y = (minY : total_Width/Ny : maxY);
    [X, Y] = meshgrid(x, y);
    
    switch selectedDirection
        case 'radiobutton_total'
            EE0 = (mag_plus(ii) * (exp(-1i*kz(ii)*(X))) + mag_minus(ii)*exp(1i*kz(ii)*(X))).*exp(-1i*ky(ii)*Y);
        case 'radiobutton_forward'
        EE0 = (mag_plus(ii) * (exp(-1i*kz(ii)*(X))) ).*exp(-1i*ky(ii)*Y);
        case 'radiobutton_backward'
         EE0 =  mag_minus(ii)*exp(1i*kz(ii)*(X)).*exp(-1i*ky(ii)*Y);
    end

    if ii == 1
        EE = [ EE0];
        
        % *exp(-1i*kz(ii)*d(ii))
    else
        EE = [EE, EE0];
    end
end

xx = (0 : total_Length/Nx : sum(thickness));

axes(handles.axes1);


imagesc(xx*1e3, y*1e3, real(EE));
colorbar;
current_axes = gca;
current_axes.YDir = 'normal';

colormap jet;
caxis([minC maxC]);

%% SHOW LABELs
set(handles.axes2, 'visible', 'on');

eps_r = eval(get(handles.edit_eps_r, 'string'));
mu_r = eval(get(handles.edit_mu_r, 'string'));
thickness = eval(get(handles.edit_thickness, 'string'))*1e-3;


nnLayer = length(eps_r);

x = eval(get(handles.edit_x, 'string'))*1e-3;
y = eval(get(handles.edit_y, 'string'))*1e-3;

minX = x(1);
maxX = x(end);

minY = y(1);
maxY = y(end);

axes(handles.axes2);
cla;
axis([minX, maxX, minY, maxY]*1e3);
colorbar off;

set(handles.axes2, 'YDir', 'Normal');

for ii = 1:nnLayer
    
    x_location = sum(thickness(1:ii));
    txt_y_location = 0;
    txt_x_location = (sum(thickness(1:ii))+sum(thickness(1:(ii-1))))/2*1e3;
    
    txt = ['Layer', num2str(ii)];
    txt_eps = ['\epsilon_r = ', num2str(eps_r(ii))];
    txt_mu = ['\mu_r = ', num2str(mu_r(ii))];
    r = 0;
    g = 0;
    b = 0;
    
    if (x_location > minX && x_location < maxX)
        
        line([x_location x_location]*1e3, [minY, maxY]*1e3, 'color', 'black', 'linestyle', '--', 'linewidth', 2);

    end
    
    text(txt_x_location, txt_y_location, txt, 'color', [r g b], 'fontsize', 10);
%     text(txt_x_location, txt_y_location-20,txt_eps, 'fontsize', 10, 'color', [r g b], 'fontweight', 'bold' );
%     text(txt_x_location, txt_y_location-40,txt_mu, 'fontsize', 10, 'color', [r g b], 'fontweight', 'bold' );
    
end
%%

% nFrame = 100;
% duration = 10;
% 
% iT = linspace(0, duration, nFrame);
% for ii = 1:nFrame
% 
% imagesc(xx, y, real(EE*exp(1i*omega*iT(ii))));
% colorbar;
% 
% colormap jet;
% caxis([minC maxC]);
% 
%  drawnow;
% 
% end