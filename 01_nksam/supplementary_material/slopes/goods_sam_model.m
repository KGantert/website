% *************************************************************************
% Shopping Time and Frictional Goods Markets:
% Implications for the New-Keynesian Business Cycle Model
% Reduced-Form NK-SaM Model Slopes
% -------------------------------------------------------------------------
% Konstantin Gantert
% Tilburg University
% Department of Economics
% k.gantert@tilburguniversity.edu
% 05/08/2025
% *************************************************************************
% -------------------------------------------------------------------------
% FULL MODEL (with and without home production)
% -------------------------------------------------------------------------
function out_sam ...
    = goods_sam_model(hhHM, cu, x, gamEH, gamSH, nuH, alpH, muM, nuM, alpM, ...
                        sig, epss, mp, ue, nuS, gamES, gamSS, pc_slope, okun, ...
                        ghh, hw_select, iHS, mp_select)

    % Steady State Calculation
    % ---------------------------------------------------------------------
    % Calculate epsW based on unemployment rate target
    epsW  = (1+ue)^nuM/((1+ue)^nuM-1);
    % FSOLVE Options
    options = optimoptions('fsolve', 'Display', 'off', ...
                            'Algorithm', 'levenberg-marquardt', ...
                            'MaxFunctionEvaluations', 20000, ...
                            'MaxIterations', 20000);
    % Steady-state for home production economy
    if hw_select == 1
        % Solve full model steady state
        steady_full_fct = ...
            @(xx) steady_full_model(xx, hhHM, cu, x, gamEH, gamSH, nuH, alpH, ...
                                    muM, nuM, alpM, sig, epss, mp, epsW, ...
                                    nuS, gamES, gamSS, ghh, mp_select);
        xx      = fsolve(steady_full_fct,[1 1],options);
        % Home production consumption
        ch    = real(xx(1));
        % Marginal utility of consumption
        uC      = real(xx(2));
        % Capacity utilization steady-state
        q       = cu;
        % Home to market consumption ratio
        chcm    = hhHM^(1-alpM)*ch^((alpM-alpH)/(1-alpH))*(1/q);
        % Market production consumption
        cm      = ch/chcm;
        % Household home production disutility parameter
        muH     = ( gamEH*(1-alpH)/ch*(gamEH*chcm^gamSH+(1-gamEH))^(1/gamSH-1) ...
                    * chcm^(gamSH-1)/(ghh+(1-ghh)/uC) )^((1-alpH)/(nuH+alpH));
        % Goods market matching efficiency
        psii    = q/((gamES*x^gamSS+(1-gamES))^(1/gamSS));
        % Composite consumption
        c    = (gamEH+(1-gamEH)*chcm^(-gamSH))^(1/gamSH) * ch;
    % Steady-state WITHOUT home production economy
    else
        steady_wo_home_fct = ...
            @(xx) steady_full_model_wo_homeprod(xx, cu, x, muM, nuM, alpM, ...
                        sig, epss, mp, epsW, nuS, gamES, gamSS, ghh, mp_select);
        xx      = fsolve(steady_wo_home_fct, [1 1], options);
        % Market production consumption
        cm      = real(xx(1));
        % Marginal utility of consumption
        uC      = real(xx(2));
        % Capacity utilization steady-state
        q       = cu;
        % Goods market tightness
        % Goods market matching efficiency
        psii    = q./((gamES.*x.^gamSS+(1-gamES)).^(1./gamSS));
        % Home production consumption
        ch      = 0;
        % Composite consumption
        c       = cm;
    end
    
    % Goods market SaM coefficients
    if mp_select == 1
        % Calculate epss based on markup target
        epss    = (q./(q-cu./mp)-gamES)./(1-gamES);
        phiGAM  = gamES.*x.^gamSS./(1-gamES);
        phiEPS  = (epss-1)./epss .* phiGAM./((1+phiGAM).*(1-iHS.*nuS));
    else
        % Calculate mp based on epss
        phiGAM  = gamES.*x.^gamSS./(1-gamES);
        phiEPS  = (epss-1)./epss .* phiGAM./((1+phiGAM).*(1-iHS.*nuS));
        mp      = phiGAM.*(1-phiEPS)./phiEPS;
    end

    % Home production coefficients
    if hw_select == 1
        chiCM   = (1-gamEH).*(cm./c).^gamSH;
        chiCH   = gamEH.*(ch./c).^gamSH;
        phiCM   = chiCM ./ (1-chiCH.*(1-gamSH-(1-ghh).*sig.*uC)./((1+nuH)./(1-alpH)-gamSH));
        phiCH   = (1-gamSH).*(1-phiCM);
    else
        chiCM   = 1;
        chiCH   = 0;
        phiCM   = chiCM;
        phiCH   = 0;
    end

    % Composite coefficients (five equation full model)
    % ---------------------------------------------------------------------

    % Utility function composite parameters
    phiUcm      = phiCM + ghh .* ( chiCM.*phiEPS./nuS.*(phiCH+(1-ghh).*sig.*phiCM) ...
                                    - chiCH./((1+nuH)./(1-alpH)-gamSH).*(1-gamSH-sig*(1-ghh)).*phiCM ...
                                    - (epsW-1)./epsW.*(1-iHS.*nuS).*phiEPS./phiGAM.*chiCM);
    phiUq       = ghh.*chiCM.*phiEPS./(psii.*phiGAM).*( 1./nuS + (epsW-1)./epsW.*(1-iHS.*nuS) );
    phiUps      = ghh.*chiCM.*phiEPS./nuS.*(1-phiEPS);
    phiUpsi     = ghh.*chiCM.*phiEPS./(psii.*phiGAM).*(1+phiGAM)./nuS;
    phiUa       = ghh.*(epsW-1)./epsW.*(1-iHS.*nuS).*phiEPS./phiGAM.*chiCM;

    % Marginal cost function composite parameters
    phiMCmc     = 1 - phiEPS.*epss./(1+phiGAM)./(1-phiEPS.*(1+(1+iHS.*nuS.*epss)./(epss-1)));
    phiMCcm     = ( (alpM+nuM)./(1-alpM)+phiCH+(1-ghh).*sig.*phiCM ) ./ phiMCmc;
    phiMCue     = nuM ./ phiMCmc ./ (1+ue);
    phiMCq      = ( (1+nuM)./(1-alpM) - phiEPS.*gamSS./phiGAM ...
                                .*(1+phiEPS./(1-phiEPS).*(1+iHS.*nuS.*epss)) ...
                                ./(1-phiEPS.*(1+(1+iHS.*nuS.*epss)./(epss-1))) ) ./ (psii.*phiMCmc);
    phiMCpsi    = ( phiEPS.*gamSS./phiGAM.*(1+phiEPS./(1-phiEPS).*(1+iHS.*nuS.*epss)) ...
                    ./(1-phiEPS.*(1+(1+iHS.*nuS.*epss)./(epss-1))) ) ./ (psii.*phiMCmc);
    phiMCa      = ((1+nuM)./(1-alpM)) ./ phiMCmc;

    % Capacity utilization composite parameters
    phiQaux     = 1 - phiEPS.*(1+(1+iHS.*nuS.*epss)./(epss-1));
    phiQq       = (1+nuS)./(psii.*phiGAM) + ((1-phiEPS).*epss./(1+phiGAM).*phiMCq)./phiQaux ...
                    - gamSS./(psii.*phiGAM).*(1-phiEPS+phiEPS.*(1+iHS.*nuS.*epss))./phiQaux;
    phiQcm      = (1-phiEPS).*epss./(1+phiGAM).*phiMCcm./phiQaux - (nuS+phiCH+(1-ghh).*sig.*phiCM);
    phiQue      = (1-phiEPS).*epss./(1+phiGAM).*phiMCue./phiQaux;
    phiQpsi     = ((1-phiEPS).*epss./(1+phiGAM)*phiMCpsi+(1-phiEPS+phiEPS.*(1+iHS.*nuS.*epss)).*gamSS./(psii.*phiGAM))./phiQaux ...
                    - (1+nuS).*(1+phiGAM)./(psii.*phiGAM);
    phiQa       = (1-phiEPS).*epss./(1+phiGAM).*phiMCa./phiQaux;

    % Capacity Utilization Slopes (Gap Variables)
    % ---------------------------------------------------------------------
    tetQcm      = phiQcm ./ phiQq;
    tetQue      = phiQue ./ phiQq;
    tetQpsi     = phiQpsi ./ phiQq;
    tetQa       = phiQa ./ phiQq;

    % Price Elasticity of Demand Slopes (Gap Variables)
    % ---------------------------------------------------------------------
    tetXIaux    = phiEPS./(1-phiEPS.*(1+(1+iHS.*nuS.*epss)./(epss-1)));
    tetXIcm     = tetXIaux .* ( epss./(1+phiGAM) .* ( phiMCcm - phiMCq.*tetQcm ) ...
                                + (1+phiEPS./(1-phiEPS).*(1+iHS.*nuS.*epss)).*gamSS./(psii.*phiGAM).*tetQcm );  
    tetXIue     = tetXIaux .* ( epss./(1+phiGAM) .* ( phiMCue - phiMCq.*tetQue ) ...
                                + (1+phiEPS./(1-phiEPS).*(1+iHS.*nuS.*epss)).*gamSS./(psii.*phiGAM).*tetQue );
    tetXIpsi    = tetXIaux .* ( epss./(1+phiGAM) .* ( phiMCpsi - phiMCq.*tetQpsi ) ...
                                + (1+phiEPS./(1-phiEPS).*(1+iHS.*nuS.*epss)).*gamSS./(psii.*phiGAM).*(1+tetQpsi) );
    tetXIa      = tetXIaux .* ( epss./(1+phiGAM) .* ( phiMCcm - phiMCq.*tetQa ) ...
                                + (1+phiEPS./(1-phiEPS).*(1+iHS.*nuS.*epss)).*gamSS./(psii.*phiGAM).*tetQa);

    % Labor share Phillips curve coefficients
    % ---------------------------------------------------------------------
    phiLSq      = 1 - phiQcm./phiQq.*phiMCq./phiMCcm;
    phiLSls     = phiQcm./phiQq.*1./phiMCcm;
    phiLSue     = phiQue./phiQq - phiQcm./phiQq.*phiMCue./phiMCcm;
    phiLSpsi    = phiQpsi./phiQq - phiQcm./phiQq.*phiMCpsi./phiMCcm;
    phiLSa      = phiQa./phiQq - phiQcm./phiQq.*phiMCa./phiMCcm;
    phiPIls     = ((1+phiGAM)./epss.*phiQaux)^(-1);
    phiPIq      = ((1-phiEPS+phiEPS.*(1+iHS.*nuS.*epss))./phiQaux-1).*gamSS./(psii.*phiGAM).*1./(1-phiEPS);

    % ---------------------------------------------------------------------
    % OUTPUT: Model Slopes and Steady-States
    % ---------------------------------------------------------------------

    % Steady-states
    % ---------------------------------------------------------------------
    out_sam.stst_ps     = phiEPS./(1-phiEPS);
    out_sam.stst_mp     = mp;
    out_sam.stst_pe    = (-epss).*(1-phiEPS);
    out_sam.stst_cm     = cm;
    out_sam.stst_mu     = chiCM.*c./cm.*(1-phiEPS).*c.^(-sig);

    % Capacity utilization slopes
    % ---------------------------------------------------------------------
    out_sam.CU_cm   = tetQcm;
    out_sam.CU_ue   = tetQue;
    out_sam.CU_psi  = tetQpsi;
    out_sam.CU_a    = tetQa;

    % Price elasticity of demand slopes
    % ---------------------------------------------------------------------
    out_sam.PE_cm   = (-1).*tetXIcm;
    out_sam.PE_ue   = (-1).*tetXIue;
    out_sam.PE_psi  = (-1).*tetXIpsi;
    out_sam.PE_a    = (-1).*tetXIa;

    % Labor wedge slopes
    % ---------------------------------------------------------------------
    out_sam.LW_cm   = phiCH + (1-ghh).*sig.*phiCM - (1-epss.*phiEPS)./(epss.*phiEPS).*tetXIcm ...
                        + gamSS./(1-phiEPS).*(1+phiGAM)./(epss.*phiGAM).*tetQcm./psii;
    out_sam.LW_ue   = nuM./(1+ue) - (1-epss.*phiEPS)./(epss.*phiEPS).*tetXIue ...
                         + gamSS./(1-phiEPS).*(1+phiGAM)./(epss.*phiGAM).*tetQue./psii;

    % Marginal costs slopes
    mc_aux          = (1-phiEPS.*((1+phiGAM)./epss.*(1-phiEPS.*(1+(1+iHS.*nuS.*epss)./(epss-1)))).^(-1)).^(-1);
    out_sam.MC_cm   = mc_aux .* ( (alpM+nuM)./(1-alpM) + phiCH + (1-ghh).*sig.*phiCM ...
                                    - ( (1+nuM)./(1-alpM) - phiEPS.*(1+phiEPS./(1-phiEPS).*(1+iHS.*nuS.*epss))./(1-phiEPS.*(1+(1+iHS.*nuS.*epss)./(epss-1))).*gamSS./phiGAM).*out_sam.CU_cm./psii);
    out_sam.MC_ue   = mc_aux .* ( nuM./(1+ue) ...
                                    - ( (1+nuM)./(1-alpM) - phiEPS.*(1+phiEPS./(1-phiEPS).*(1+iHS.*nuS.*epss))./(1-phiEPS.*(1+(1+iHS.*nuS.*epss)./(epss-1))).*gamSS./phiGAM).*out_sam.CU_ue./psii);

    % Euler equation slope
    % ---------------------------------------------------------------------
    out_sam.AD_cm   = phiCH + tetXIcm ...
                        + sig/uC .* ( phiUcm + phiUq.*tetQcm - phiUps./phiEPS.*tetXIcm );
    out_sam.AD_ue   = tetXIue + sig/uC .* ( phiUq.*tetQue - phiUps./phiEPS.*tetXIue );

    % Phillips curve slope
    % ---------------------------------------------------------------------
    out_sam.kapP    = 1/pc_slope .* phiGAM.*(1-phiEPS)./(phiEPS.*(1-iHS.*nuS)).*(epss-1)./(epss) ...
                        .* ( phiPIls + phiPIq.*phiLSls./phiLSq );
    out_sam.AS_cm   = (1+phiGAM)./out_sam.kapP .* ( (1-phiEPS)./(phiEPS.*(1-iHS.*nuS)).*tetXIcm - 1./(1-iHS.*nuS).*gamSS./phiGAM.*tetQcm );
    out_sam.AS_ue   = (1+phiGAM)./out_sam.kapP .* ( (1-phiEPS)./(phiEPS.*(1-iHS.*nuS)).*tetXIue - 1./(1-iHS.*nuS).*gamSS./phiGAM.*tetQue );

    % Real wage qquation coefficients
    % ---------------------------------------------------------------------
    out_sam.WG_cm   = phiMCcm - alpM./(1-alpM) + (1./(1-alpM)-phiMCq).*tetQcm;
    out_sam.WG_ue   = phiMCue + (1/(1-alpM)-phiMCq)*tetQue;

    % Applying Okun's Law
    % ---------------------------------------------------------------------
    out_sam.CU_ag   = out_sam.CU_cm + out_sam.CU_ue./okun;
    out_sam.PE_ag   = out_sam.PE_cm + out_sam.PE_ue./okun;
    out_sam.LW_ag   = out_sam.LW_cm + out_sam.LW_ue./okun;
    out_sam.AS_ag   = out_sam.AS_cm + out_sam.AS_ue./okun;
    out_sam.AD_ag   = out_sam.AD_cm + out_sam.AD_ue./okun;
    out_sam.WG_ag   = out_sam.WG_cm + out_sam.WG_ue./okun;
    out_sam.MC_ag   = out_sam.MC_cm + out_sam.MC_ue./okun;

    % ---------------------------------------------------------------------
    % Steady state function full model
    function F = steady_full_model(xx, hhHM, cu, x, gamEH, gamSH, nuH, alpH, ...
                                    muM, nuM, alpM, sig, epss, mp, epsW, ...
                                    nuS, gamES, gamSS, ghh, mp_select)
        % ch = xx(1), uC = xx(2)
        % Capacity utilization steady-state
        q       = cu;
        % Home to market consumption ratio
        chcm    = hhHM^(1-alpM)*xx(1)^((alpM-alpH)/(1-alpH))*(1/q);
        % Household home production disutility parameter
        muH   = ( gamEH*(1-alpH)/xx(1)*(gamEH*chcm^gamSH+(1-gamEH))^(1/gamSH-1) ...
                  * chcm^(gamSH-1)*1/(ghh+(1-ghh)/xx(2)) )^((1-alpH)/(nuH+alpH));
        % Goods market matching efficiency
        psii    = q/((gamES*x^gamSS+(1-gamES))^(1/gamSS));
        % Composite consumption
        c    = (gamEH+(1-gamEH)*chcm^(-gamSH))^(1/gamSH) * xx(1);
        % Goods market SaM coefficients
        if mp_select == 1
            % Calculate epss based on markup target
            epss    = (q/(q-cu/mp)-gamES)/(1-gamES);
            phiGAM  = gamES.*x.^gamSS./(1-gamES);
            phiEPS  = (epss-1)./epss * phiGAM./((1+phiGAM).*(1-iHS.*nuS));
        else
            % Calculate mp based on epss
            phiGAM  = gamES.*x.^gamSS./(1-gamES);
            phiEPS  = (epss-1)./epss * phiGAM./((1+phiGAM).*(1-iHS.*nuS));
            mp      = phiGAM.*(1-phiEPS)./phiEPS;
        end
        % Auxiliary equation for uC
        csCM = (epss-1)/epss*(gamES*x^gamSS)/(gamES*x^gamSS+(1-gamES))*(1-gamEH)*c^(1-gamSH)*xx(1)^gamSH*chcm^(-gamSH);
        % Solve equilibrium conditions for cm and ch
        F(1) = (xx(1)/chcm)^((1+nuM)/(1-alpM)-1) * (ghh+(1-ghh)/xx(2)) ...
                - (psii*(gamES*x^gamSS+(1-gamES))^(1/gamSS))^((1+nuM)/(1-alpM))*(1-alpM)/muM * (epsW-1)/epsW * (epss-1)/(epss*(1+gamES*x^gamSS/(1-gamES)))*(1-gamEH)*(gamEH*chcm^gamSH+(1-gamEH))^(1/gamSH-1);
        F(2) = xx(2)^(-1/sig) - ( c - ghh/(1+nuS)*csCM ...
                         - ghh*muH/(1+nuH)*xx(1)^((1+nuH)/(1-alpH)) ...
                         - ghh*muM/(1+nuM)*(xx(1)/(psii*chcm))^((1+nuM)/(1-alpM)) );
    end
    % ---------------------------------------------------------------------
    % Steady state function full model without homeproduction
    function F = steady_full_model_wo_homeprod(xx, cu, x, ...
                                    muM, nuM, alpM, sig, epss, mp, epsW, ...
                                    nuS, gamES, gamSS, ghh, mp_select)
        % cm = xx(1), uC = xx(2)
        % Capacity utilization steady-state
        q       = cu;
        % Goods market matching efficiency
        psii    = q/((gamES*x^gamSS+(1-gamES))^(1/gamSS));
        % Goods market SaM coefficients
        if mp_select == 1
            % Calculate epss based on markup target
            epss    = (q/(q-cu/mp)-gamES)/(1-gamES);
            phiGAM  = gamES.*x.^gamSS./(1-gamES);
            phiEPS  = (epss-1)./epss * phiGAM./((1+phiGAM).*(1-iHS.*nuS));
        else
            % Calculate mp based on epss
            phiGAM  = gamES.*x.^gamSS./(1-gamES);
            phiEPS  = (epss-1)./epss * phiGAM./((1+phiGAM).*(1-iHS.*nuS));
            mp      = phiGAM.*(1-phiEPS)./phiEPS;
        end
        % Auxiliary equation for uC
        csCM = (epss-1)/epss*(gamES*x^gamSS)/(gamES*x^gamSS+(1-gamES))*xx(1);
        % Solve equilibrium conditions for cm and ch
        F(1) = xx(1)^((1+nuM)/(1-alpM)-1) * (ghh+(1-ghh)/xx(2)) ...
                - (psii*(gamES*x^gamSS+(1-gamES))^(1/gamSS))^((1+nuM)/(1-alpM))*(1-alpM)/muM ...
                    * (epsW-1)/epsW * (epss-1)/(epss*(1+gamES*x^gamSS/(1-gamES)));
        F(2) = xx(2)^(-1/sig) - ( xx(1) - ghh/(1+nuS)*csCM ...
                         - ghh*muM/(1+nuM)*(xx(1)/psii)^((1+nuM)/(1-alpM)) );
    end

end