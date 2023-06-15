function [eps_r, thickness] = eval_eps_r_thickness_from_mu_r(handles)

mu_r = eval(get(handles.edit_mu_r, 'string'));
eps_r0 = eval(get(handles.edit_eps_r, 'string'));
thickness0 = eval(get(handles.edit_thickness, 'string'));

if (length(mu_r) > length(eps_r0))
    n = length(mu_r) - length(eps_r0);
    eps_r = [eps_r0; ones(n, 1)];
    
    %%
    if (length(mu_r) > length(thickness0))
        n = length(mu_r) - length(thickness0);
        thickness = [thickness0; ones(n, 1)];
        
    elseif  (length(mu_r) < length(thickness0))
        n = length(thickness0) - length(mu_r);
        thickness0((end - n+1):end) = [];
        thickness = thickness0;
    else
        thickness = thickness0;
    end
    
    %%
elseif (length(mu_r) < length(eps_r0))
    n = length(eps_r0) - length(mu_r);
    eps_r0((end - n+1):end) = [];
    eps_r = eps_r0;
    
    %%
    if (length(mu_r) > length(thickness0))
        n = length(mu_r) - length(thickness0);
        thickness = [thickness0; ones(n, 1)];
        
    elseif  (length(mu_r) < length(thickness0))
        n = length(thickness0) - length(mu_r);
        thickness0((end - n+1):end) = [];
        thickness = thickness0;
    else
        thickness = thickness0;
    end
    
    %%
    
else
    eps_r = eps_r0;
    
    %%
    if (length(mu_r) > length(thickness0))
        n = length(mu_r) - length(thickness0);
        thickness = [thickness0; ones(n, 1)];
        
    elseif  (length(mu_r) < length(thickness0))
        n = length(thickness0) - length(mu_r);
        thickness0((end - n+1):end) = [];
        thickness = thickness0;
    else
        thickness = thickness0;
    end
    
    %%
    
end



end