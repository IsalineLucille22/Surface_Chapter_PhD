function [Evol_Col, Biomass, mu_eff] = Track_Evol_Col(num_col, Mass_Cell, mu_evol, Res_Cons_ind)
%Return the number of species, the biomass and the mean of the maximum
%growth rate for the resource "Res_Cons_Ind) for each mirco-colony at a time t_i for the species Species_index
nb_col = unique(num_col);
N_col = max(nb_col);
[Evol_Col, Biomass, mu_eff] = deal(zeros(1, N_col));
for i=1:N_col
    Evol_Col(i) = sum(num_col == i);
    Biomass(i) = sum(Mass_Cell(num_col == i));
    mu_eff(i) = mean(mu_evol(Res_Cons_ind, num_col == i));
end
