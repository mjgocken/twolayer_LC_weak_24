# twolayer_LC_weak_24
Two layer model of the tear film. LL is represented as LC with weak elasticity; AL is Newtonian.   
These files were used for the La Matematica submission; figures mentioned refer to those in the paper.

These Matlab files model the tear film in two layers (LL is LC; AL is Newtonian).
   - [finite difference code](twolayer_withLC_FD_mjg.m)
   - Fig 4: extrema v theta_B
        * [code](extrema_v_thetaB.m)
        * [data](LC_3times5CB_thetaDepEvap_extrema_v_theta_dim.mat)
   - Fig 5 and 6: final time and extrema through time
        * [code](comparison_plots.m)
        * [data: 0](LC_3times5CB_thetaDepEvap_thetaB0_dim_finer.mat)
        * [data: pi/18](LC_3times5CB_thetaDepEvap_thetaBpi18_dim_finer.mat)
        * [data: pi/12](LC_3times5CB_thetaDepEvap_thetaBpi12_dim_finer.mat)
        * [data: pi/4](LC_3times5CB_thetaDepEvap_thetaBpi4_dim_finer.mat)
        * [data: pi/2](LC_3times5CB_thetaDepEvap_thetaBpi2_dim_finer.mat)<br>
The Mathematica scripts contain the equations at each order of epsilon using the asymptotic expansions of the dependent variables.<br>
    - [expansions](twolayer_expansions.nb)  
    - [order one](twolayer_order1.nb)  
    - [order epsilon](twolayer_orderepsilon.nb)  
    - [order epsilon squared](twolayer_orderepsilonsquared.nb)  
    - [osmolarity equations](salt.nb)  
    - [simplificiations](simplify_Fprime.nb)  
    - [u_21](for_u21_feb7.nb)  
    - [u_20](for_u20_using_u22yy_subforu21y_feb7.nb)  
