function pi_name_red = reduce_pi_name(pi_name);

pi_name_red='';
[u,r]=strtok(pi_name,' ');
while isempty(r)==0

    u_upper=upper(u);
    u_lower=lower(u);
    if (strcmp(u,u_upper)==0)&(strcmp(u,u_lower)==0); % mot en minuscule, avec capitale, on ne garde que l'initiale
        pi_name_red=[pi_name_red, u(u==u_upper),'. '];
    else (strcmp(u,u_lower)==1)| (strcmp(u,u_upper)==1);  % mot en minuscupe ou majuscule, on garde
        if (strcmp(u,u_upper)==1)
        u=u_lower;
        u(1)=u_upper(1);
        end
        pi_name_red=[pi_name_red, u ,' '];
    
    end
    [u,r]=strtok(r,' ');
end
u_upper=upper(u);
u_lower=lower(u);


if (strcmp(u,u_upper)==0) & (strcmp(u,u_lower)==0); % mot en minuscule, avec capitale, on ne garde que l'initiale
    pi_name_red=[pi_name_red, u(u==u_upper),'. '];
    
else (strcmp(u,u_lower)==1)|( strcmp(u,u_upper)==1);  % mot en minuscupe ou majuscule, on garde
 
   if (strcmp(u,u_upper)==1);
       u=u_lower;
       u(1)=u_upper(1);
    end
   pi_name_red=[pi_name_red, u ,' '];
end

pi_name_red = strrep(pi_name_red,' et ','/');
pi_name_red = strrep(pi_name_red,'. ','.');