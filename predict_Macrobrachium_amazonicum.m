%% predict_Macrobrachium_amazonicum
% Obtains predictions, using parameters and data

%%%%%%NOTE: This script refers only to females of the species, need to adapt it for male morphotypes%%%%%%

%%
function [prdData, info] = predict_Macrobrachium_amazonicum_female(par, data, auxData)
  % created by Starrlight Augustine, Dina Lika, Bas Kooijman, Goncalo Marques and Laure Pecquerie 2015/01/30; 
  % last modified 2015/07/29
  
  %% Syntax
  % [prdData, info] = <../predict_my_pet.m *predict_my_pet*>(par, data, auxData)
  
  %% Description
  % Obtains predictions, using parameters and data
  %
  % Input
  %
  % * par: structure with parameters (see below)
  % * data: structure with data (not all elements are used)
  % * auxData : structure with temp data and other potential environmental data
  %  
  % Output
  %
  % * prdData: structure with predicted values for data
  % * info: identified for correct setting of predictions (see remarks)
  
  %% Remarks
  % Template for use in add_my_pet.
  % The code calls <parscomp_st.html *parscomp_st*> in order to compute
  % scaled quantities, compound parameters, molecular weights and compose
  % matrixes of mass to energy couplers and chemical indices.
  % With the use of filters, setting info = 0, prdData = {}, return, has the effect
  % that the parameter-combination is not selected for finding the
  % best-fitting combination; this setting acts as customized filter.
  
  %% Example of a customized filter
  % See the lines just below unpacking
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
    
  % compute temperature correction factors
  TC= tempcorr(temp.ab, T_ref, T_A); 
  TC_WwN = tempcorr(temp.WwN, T_ref, T_A);
  TC_LN = tempcorr(temp.LN, T_ref, T_A);
  kT_M = k_M * TC;                  % 1/d, som maint rate coeff
 
% uncomment if you need this for computing moles of a gas to a volume of gas
% - else feel free to delete  these lines
% molar volume of gas at 1 bar and 20 C is 24.4 L/mol
% T = C2K(20); % K, temp of measurement equipment- apperently this is
% always the standard unless explicitely stated otherwise in a paper (pers.
% comm. Mike Kearney).
% X_gas = T_ref/ T/ 24.4;  % M,  mol of gas per litre at T_ref and 1 bar;
  
% zero-variate data

  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

  % birth
   L_b = L_m * l_b;                  % cm, structural length at birth at f
   Ww_b = L_b^3 *(1 + f * w);        % g, wet weight at birth
   Wd_b = L_b^3 * d_V * (1 + f * w); % g, dry weight at birth
   aT_b = t_0 + t_b/ k_M/ TC;        % d, age at birth
  
  % metam
  L_j = L_m * l_j;                  % cm, structural length at metam at f
  Lw_j = L_j/ del_M;                % cm, total length at metam
  aT_j = (t_j-t_b)/ kT_M;           % d, time since birth at metam at f and T

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, total length at puberty at f
  Ww_p = L_p^3 *(1 + f * w);        % g, wet weight at puberty 
  aT_p = (t_p-t_j)/ kT_M;           % d, time since metam at puberty at f
  

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
  RT_i = TC * reprod_rate_j(L_i, f, pars_R);                    % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC;               % d, mean life span at T
  
  % pack to output
  % the names of the fields in the structure must be the same as the data names in the mydata file
  prdData.ab = aT_b;
  prdData.tj = aT_j;
  prdData.ap = aT_p;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wdb = Wd_b;
  prdData.Ri = RT_i;
  
  % uni-variate data
  
  % wet weight-fecundity
  L = (WwN(:,1)./ (1 + f * w)).^(1/3);   % cm, structural length
  WwE = TC_WwN * 365 * reprod_rate_j(L, f, pars_R);  % annual fecundity
  
  % length-wet weight
  EWw = (LW(:,1) * del_M).^3 * (1 + f * w); % g, wet weight  
  
  % length-fecundity
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
  EN = TC_LN * 365 * reprod_rate_j(LN(:,1) * del_M, f, pars_R);% #, fecundity and length
  
  % pack to output
  % the names of the fields in the structure must be the same as the data names in the mydata file
  prdData.WwN = WwE;
  prdData.LW = EWw;
  prdData.LN = EN;  
end 
