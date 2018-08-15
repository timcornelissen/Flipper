function U = solve_induced(Grids,In,indices,grid,ind,filter)
% SOLVE_INDUCED solves the self-consistent induced dipole problem
% [U] = SOLVE_INDUCED(Grids,In,indices,grid,ind,filter)
% indices are the indices of the dipole to include in the calculation, grid contains their orientation
% ind and filter are filters to prevent taking energy terms outside of the interaction range into account

dipoles = bsxfun(@times,grid(Grids.neighbourindex(indices,:)),Grids.dipoles(indices,:,:));
dipolestotal = dipoles;

for k=1:In.iter
    field = -bsxfun(@times,Grids.r(indices,:),dipolestotal-3*bsxfun(@times,dot(dipolestotal,Grids.ru(indices,:,:),3),Grids.ru(indices,:,:)));
    reactionfield = In.reaction*min(Grids.r(1,end))*dipolestotal; %Take maximum range of current dipole
    induced_single = In.pol*In.ke*(squeeze(sum(field+reactionfield,2))+ [0 0 In.Eappl/In.ke/In.mu]); %Induced dipole
    dipolestotal = dipoles+induced_single(ind).*repmat(filter,1,1,3); %remove dipoles outside range
end

dip = grid(indices).*Grids.posgrid(indices,:,2)+induced_single;
U = In.fieldconstant*sum(dot(dip,squeeze(sum(field+reactionfield,2))+ [0 0 In.Eappl/In.ke/In.mu],2));
end