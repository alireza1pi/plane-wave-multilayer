function build_structure(nLayer, eps_r, mu_r, thickness, lastLayerFlag, mws )

thickness = thickness(1:nLayer);
eps_r = eps_r(1:nLayer);
mu_r = mu_r(1:nLayer);

xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;

%Create The ii'th Layer

for ii = 1:nLayer
        
%     d = sum(thickness(1:ii)); %
    if ii ==1
        zmin = 0;
    else
        zmin = 1000*sum(thickness(1:ii-1));
    end
    
    zmax = 1000*sum(thickness(1:ii));
    
   
    
    %Define Material
    material_Name = strcat('Material_Layer', num2str(ii) );
    
    define_material(material_Name, eps_r(ii), mu_r(ii), 0, 0, mws);
    
    brick = mws.invoke('Brick');
    brick.invoke('Reset');
    
    brick_Name = strcat('Brick_Layer', num2str(ii));
    
    brick.invoke('Name', brick_Name);
    brick.invoke('Component', 'component1');
    brick.invoke('Material', material_Name);
    brick.invoke('Xrange', num2str(xmin), num2str(xmax));
    brick.invoke('Yrange', num2str(ymin), num2str(ymax));
    brick.invoke('Zrange', num2str(-zmin), num2str(-zmax));
    
    brick.invoke('Create');
    
    
    
    
    

end