function [mu_r, eps_r] = eval_mu_r_eps_r_from_thickness(handles)

thickness = eval(get(handles.edit_thickness, 'string'));
eps_r0 = eval(get(handles.edit_eps_r, 'string'));
mu_r0 = eval(get(handles.edit_mu_r, 'string'));

if (length(thickness) > length(eps_r0))
    n = length(thickness) - length(eps_r0);
    eps_r = [eps_r0; ones(n, 1)];
    
    %%
    if (length(thickness) > length(mu_r0))
        n = length(thickness) - length(mu_r0);
        mu_r = [mu_r0; ones(n, 1)];
        
    elseif  (length(thickness) < length(mu_r0))
        n = length(mu_r0) - length(thickness);
        mu_r0((end - n+1):end) = [];
        mu_r = mu_r0;
    else
        mu_r = mu_r0;
    end
    
    %%
elseif (length(thickness) < length(eps_r0))
    n = length(eps_r0) - length(thickness);
    eps_r0((end - n+1):end) = [];
    eps_r = eps_r0;
    
    %%
    if (length(thickness) > length(mu_r0))
        n = length(thickness) - length(mu_r0);
        mu_r = [mu_r0; ones(n, 1)];
        
    elseif  (length(thickness) < length(mu_r0))
        n = length(mu_r0) - length(thickness);
        mu_r0((end - n+1):end) = [];
        mu_r = mu_r0;
    else
        mu_r = mu_r0;
    end
    
    %%
    
else
    eps_r = eps_r0;
    
    %%
    if (length(thickness) > length(mu_r0))
        n = length(thickness) - length(mu_r0);
        mu_r = [mu_r0; ones(n, 1)];
        
    elseif  (length(thickness) < length(mu_r0))
        n = length(mu_r0) - length(thickness);
        mu_r0((end - n+1):end) = [];
        mu_r = mu_r0;
    else
        mu_r = mu_r0;
    end
    
    %%
    
end



end