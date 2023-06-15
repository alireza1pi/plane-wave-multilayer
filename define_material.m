function define_material(name, eps_r, mu_r, sigma, tanD, mws)

R = rand(1);
G = rand(1);
B = rand(1);
material = mws.invoke('Material');
material.invoke('Reset');
material.invoke('Name', name);
material.invoke('Folder', '');
material.invoke('Type', 'Normal');
material.invoke('Epsilon', num2str(eps_r));
material.invoke('Mu', num2str(mu_r));
material.invoke('Sigma', num2str(sigma));
material.invoke('TanD', num2str(tanD));
material.invoke('TanDModel', 'ConstTanD');
material.invoke('Colour', num2str(R), num2str(G), num2str(B));

material.invoke('Create');

end