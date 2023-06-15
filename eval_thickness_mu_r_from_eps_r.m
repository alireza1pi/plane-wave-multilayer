function [thickness, mu_r] = eval_thickness_mu_r_from_eps_r(handles)

eps_r = eval(get(handles.edit_eps_r, 'string'));
thickness0 = eval(get(handles.edit_thickness, 'string'));
mu_r0 = eval(get(handles.edit_mu_r, 'string'));

if (length(eps_r) > length(thickness0))
    n = length(eps_r) - length(thickness0);
    thickness = [thickness0; zeros(n, 1)];
    
    %%
    if (length(eps_r) > length(mu_r0))
        n = length(eps_r) - length(mu_r0);
        mu_r = [mu_r0; ones(n, 1)];
        
    elseif  (length(eps_r) < length(mu_r0))
        n = length(mu_r0) - length(eps_r);
        mu_r0((end - n+1):end) = [];
        mu_r = mu_r0;
    else
        mu_r = mu_r0;
    end
    
    %%
elseif (length(eps_r) < length(thickness0))
    n = length(thickness0) - length(eps_r);
    thickness0((end - n+1):end) = [];
    thickness = thickness0;
    
    %%
    if (length(eps_r) > length(mu_r0))
        n = length(eps_r) - length(mu_r0);
        mu_r = [mu_r0; ones(n, 1)];
        
    elseif  (length(eps_r) < length(mu_r0))
        n = length(mu_r0) - length(eps_r);
        mu_r0((end - n+1):end) = [];
        mu_r = mu_r0;
    else
        mu_r = mu_r0;
    end
    
    %%
    
else
    thickness = thickness0;
    
    %%
    if (length(eps_r) > length(mu_r0))
        n = length(eps_r) - length(mu_r0);
        mu_r = [mu_r0; ones(n, 1)];
        
    elseif  (length(eps_r) < length(mu_r0))
        n = length(mu_r0) - length(eps_r);
        mu_r0((end - n+1):end) = [];
        mu_r = mu_r0;
    else
        mu_r = mu_r0;
    end
    
    %%
    
end



end