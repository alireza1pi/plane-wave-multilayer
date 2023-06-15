function animate_field(EE, thickness, nFrame, duration, handles)

freq = str2double(get(handles.edit_freq, 'string'))*1e9;
omega = 2*pi*freq;

nnLayer = size(thickness, 2);

thickness = thickness(1:nnLayer);
total_Length = sum(thickness);

x = eval(get(handles.edit_x, 'string'));
y = eval(get(handles.edit_x, 'string'));
% Nx = size(x, 2);

minC = str2double(get(handles.edit_minC, 'string'));
maxC = str2double(get(handles.edit_maxC, 'string'));

% x = 0 : total_Length/Nx : sum(thickness);

% y = linspace(minY, maxY, nY);

axes(handles.axes1);

iT = linspace(0, duration, nFrame);
for ii = 1:nFrame

imagesc(x, y, real(EE*exp(1i*omega*iT(ii))));
set(handles.axes1, 'YDir', 'normal');
colorbar;

colormap jet;
caxis([minC maxC]);

 drawnow;

end

end