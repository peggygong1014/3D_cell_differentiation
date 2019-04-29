function [Perm_tissue,P_tissue,E_tissue]=AutoMP_Select(DiffTypeN)
global E_granulation E_fibrous E_cartilage E_marrow E_immature E_mature E_cortical
global Perm_granulation Perm_fibrous Perm_cartilage Perm_marrow Perm_immature Perm_mature Perm_cortical
global P_granulation P_mature P_cortical P_fibrous P_cartilage P_marrow P_immature
switch DiffTypeN
    case 0 %Resorption
        E_tissue=0;
        Perm_tissue = 0;
        P_tissue = 0;
    case 1 %Mature
        E_tissue=E_mature;
        Perm_tissue = Perm_mature;
        P_tissue = P_mature;
    case 2
        E_tissue=E_immature;
        Perm_tissue = Perm_immature;
        P_tissue = P_immature;
    case 3
        E_tissue=E_cartilage;
        Perm_tissue = Perm_cartilage;
        P_tissue = P_cartilage;
    case 4
        E_tissue=E_fibrous;
        Perm_tissue = Perm_fibrous;
        P_tissue = P_fibrous;
    case 5
        E_tissue=E_granulation;
        Perm_tissue = Perm_granulation;
        P_tissue = P_granulation;
end