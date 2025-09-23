% *************************************************************************
% Shopping Time and Frictional Goods Markets:
% Implications for the New-Keynesian Business Cycle Model
% Reduced-Form Baseline NK Model Slopes
% -------------------------------------------------------------------------
% Konstantin Gantert
% Tilburg University
% Department of Economics
% k.gantert@tilburguniversity.edu
% 05/08/2025
% *************************************************************************
% ---------------------------------------------------------------------
% BASELINE MODEL (without goods market SaM)
% ---------------------------------------------------------------------
function out_nk ...
    = baseline_model(hhHM, gamEH, gamSH, nuH, alpH, muM, nuM, alpM, sig, ...
                        epss, mp, ue, pc_slope, okun, ghh, hw_select, mp_select)

    % Steady State Calculation
    % ---------------------------------------------------------------------
    if mp_select == 1
        % Calculate epss based on markup target
        epss    = 1/(1-1/mp);
    else
        % Calculate mp based on epss
        mp = epss/(epss-1);
    end
    % Calculate epsW based on unemployment rate target
    epsW  = (1+ue)^nuM/((1+ue)^nuM-1);
    % FSOLVE Options
    options = optimoptions('fsolve', 'Display', 'off', ...
                            'Algorithm', 'levenberg-marquardt', ...
                            'MaxFunctionEvaluations', 20000, ...
                            'MaxIterations', 20000);
    % Steady state calculation with home production
    if hw_select == 1
        % Solve full model steady state
        steady_nk_fct = ...
            @(xx) steady_nk_model(xx, hhHM, gamEH, gamSH, nuH, alpH, ...
                                    muM, nuM, alpM, sig, epss, epsW, ghh);
        xx      = fsolve(steady_nk_fct,[1 1],options);
        % Home production consumption
        ch    = real(xx(1));
        % Marginal utility of consumption
        uC    = real(xx(2));
        % Home to market consumption ratio
        chcm  = hhHM^(1-alpM)*ch^((alpM-alpH)/(1-alpH));
        % Household home production disutility parameter
        muH     = ( gamEH*(1-alpH)/ch*(gamEH*chcm^gamSH+(1-gamEH))^(1/gamSH-1) ...
                    * chcm^(gamSH-1)*1/(ghh+(1-ghh)/uC) )^((1-alpH)/(nuH+alpH));
        % Market production consumption
        cm    = ch/chcm;
        % Composite consumption
        c     = (gamEH*ch^gamSH+(1-gamEH)*cm^gamSH)^(1/gamSH);
    % Steady state calculation WITHOUT home production
    else
        % Marginal utility steady state
        uC    = (1-ghh) * ( (1-alpM)/(muM)*(epsW-1)/epsW*(epss-1)/epss )^(((1-alpM)/(nuM+alpM))/((-1)/sig-(1-alpM)/(nuM+alpM))) ...
                    + ghh * ( (1-(1-alpM)/(1+nuM)*(epsW-1)/epsW*(epss-1)/epss) ...
                            * ((1-alpM)/(muM)*(epsW-1)/epsW*(epss-1)/epss)^((1-alpM)/(nuM+alpM)) )^(-sig);
        % Market consumption steady state (NK model)
        cm    = (1-ghh) * ((1-alpM)/(muM)*(epsW-1)/epsW*(epss-1)/epss*uC)^((1-alpM)/(nuM+alpM)) ...
                    + ghh * ((1-alpM)/(muM)*(epsW-1)/epsW*(epss-1)/epss)^((1-alpM)/(nuM+alpM));
        % Composite consumption steady state (NK model)
        c     = cm;
    end
    
    % Home production coefficients
    if hw_select == 1
        chiCM   = (1-gamEH).*(cm./c).^gamSH;
        chiCH   = gamEH.*(ch./c).^gamSH;
        phiCM   = chiCM ./ (1-chiCH.*(1-gamSH-(1-ghh).*sig.*uC)./((1+nuH)./(1-alpH)-gamSH));
        phiCH   = (1-gamSH)*(1-phiCM);
    else
        chiCM   = 1;
        chiCH   = 0;
        phiCM   = 1;
        phiCH   = 0;
    end

    % Composite coefficients (three equation full model)
    % ---------------------------------------------------------------------

    % Utility function composite parameters
    phiUcm      = phiCM + ghh .* ( - chiCH./((1+nuH)./(1-alpH)-gamSH).*(1-gamSH-sig*(1-ghh)).*phiCM ...
                                    - (epsW-1)./epsW.*(epss-1)./epss.*chiCM);
    phiUa       = ghh.*(epsW-1)./epsW.*(epss-1)./epss.*chiCM;

    % Marginal cost function composite parameters
    phiMCmc     = 1;
    phiMCcm     = ( (alpM+nuM)./(1-alpM)+phiCH+(1-ghh).*sig.*phiCM ) ./ phiMCmc;
    phiMCue     = nuM ./ phiMCmc ./ (1+ue);
    phiMCq      = ( (1+nuM)./(1-alpM) ) ./ phiMCmc;
    phiMCa      = ((1+nuM)./(1-alpM)) ./ phiMCmc;

    % ---------------------------------------------------------------------
    % Output: Model Slopes and Steady-States
    % ---------------------------------------------------------------------

    % Steady-states
    % ---------------------------------------------------------------------
    out_nk.stst_ps    = 0;
    out_nk.stst_mp    = mp;
    out_nk.stst_pe    = (-epss);
    out_nk.stst_cm    = cm;
    out_nk.stst_mu    = chiCM.*c./cm.*c.^(-sig);

    % Labor wedge slopes
    % ---------------------------------------------------------------------
    out_nk.LW_cm    = phiCH + (1-ghh).*sig.*phiCM - ( phiMCcm );
    out_nk.LW_ue    = nuM./(1+ue) - ( phiMCue );

    % Marginal costs slopes
    % ---------------------------------------------------------------------
    out_nk.MC_cm    = ((alpM+nuM)./(1-alpM)+phiCH+(1-ghh).*sig.*phiCM);
    out_nk.MC_ue    = nuM./(1+ue);

    % Euler equation slope
    % ---------------------------------------------------------------------  
    out_nk.AD_cm    = phiCH + sig/uC .* ( phiUcm );
    out_nk.AD_ue    = 0;

    % Phillips curve slope
    % ---------------------------------------------------------------------
    out_nk.kapP     = epss/pc_slope;
    out_nk.AS_cm    = 1./out_nk.kapP .* ( epss .* phiMCcm );
    out_nk.AS_ue    = 1./out_nk.kapP .* ( epss .* phiMCue );

    % Real wage equation coefficients
    % ---------------------------------------------------------------------
    out_nk.WG_cm    = phiMCcm - alpM./(1-alpM);
    out_nk.WG_ue    = phiMCue;

    % Applying Okun's Law
    % ---------------------------------------------------------------------
    out_nk.AS_ag    = out_nk.AS_cm + out_nk.AS_ue./okun;
    out_nk.AD_ag    = out_nk.AD_cm + out_nk.AD_ue./okun;
    out_nk.WG_ag    = out_nk.WG_cm + out_nk.WG_ue./okun;
    out_nk.LW_ag    = out_nk.LW_cm + out_nk.LW_ue./okun;
    out_nk.MC_ag    = out_nk.MC_cm + out_nk.MC_ue./okun;


    % ---------------------------------------------------------------------
    % Steady state function full model
    function F = steady_nk_model(xx, hhHM, gamEH, gamSH, nuH, alpH, ...
                                    muM, nuM, alpM, sig, epss, epsW, ghh)
        % ch/cm = xx(1), uC = xx(2)
        % Homeproduction consumption
        chcm    = hhHM^(1-alpM)*xx(1)^((alpM-alpH)/(1-alpH));
        % Household search disutility parameter
        muH   = ( gamEH*(1-alpH)/xx(1)*(gamEH*chcm^gamSH+(1-gamEH))^(1/gamSH-1) ...
                    * chcm^(gamSH-1)*1/(ghh+(1-ghh)/xx(2)) )^((1-alpH)/(nuH+alpH));
        % Composite consumption
        c    = (gamEH+(1-gamEH)*chcm^(-gamSH))^(1/gamSH) * xx(1);
        % Solve equilibrium conditions for cm and ch
        F(1) = (xx(1)/chcm)^((1+nuM)/(1-alpM)-1) * (ghh+(1-ghh)/xx(2)) ...
                - (1-alpM)/muM * (epsW-1)/epsW*(epss-1)/epss*(1-gamEH)*(gamEH*chcm^gamSH+(1-gamEH))^(1/gamSH-1);
        F(2) = xx(2)^(-1/sig) - ( c - ghh*muH/(1+nuH)*xx(1)^((1+nuH)/(1-alpH)) ...
                         - ghh*muM/(1+nuM)*(xx(1)/(chcm))^((1+nuM)/(1-alpM)) );
    end

end