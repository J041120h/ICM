function dydt = VEGFAbeqns(t,y,p,n)
% Equations: Simplified model for Stefanini et al
% one isoform of VEGF; one receptor; one antibody
% three tissues (blood (central), tumor, rest of body)
% 2025-11-18

% t, y : current time and concentrations of all components
% p : parameter array
% n : index array (for ease of reading equations)

dydt = zeros(15,1);

% VEGFR2 and VEGF-VEGFR2 complex
dydt(n.VEGFR2_b) = p.VEGFR2prod_b - p.kintR2*y(n.VEGFR2_b) - p.konVR*y(n.VEGF_b)*y(n.VEGFR2_b) + p.koffVR*y(n.VEGFVEGFR2_b);
dydt(n.VEGFR2_t) = p.VEGFR2prod_t - p.kintR2*y(n.VEGFR2_t) - p.konVR*y(n.VEGF_t)*y(n.VEGFR2_t) + p.koffVR*y(n.VEGFVEGFR2_t);
dydt(n.VEGFR2_r) = p.VEGFR2prod_r - p.kintR2*y(n.VEGFR2_r) - p.konVR*y(n.VEGF_r)*y(n.VEGFR2_r) + p.koffVR*y(n.VEGFVEGFR2_r);

dydt(n.VEGFVEGFR2_b) =       - p.kintVR2*y(n.VEGFVEGFR2_b) + p.konVR*y(n.VEGF_b)*y(n.VEGFR2_b) - p.koffVR*y(n.VEGFVEGFR2_b);
dydt(n.VEGFVEGFR2_t) =       - p.kintVR2*y(n.VEGFVEGFR2_t) + p.konVR*y(n.VEGF_t)*y(n.VEGFR2_t) - p.koffVR*y(n.VEGFVEGFR2_t);
dydt(n.VEGFVEGFR2_r) =       - p.kintVR2*y(n.VEGFVEGFR2_r) + p.konVR*y(n.VEGF_r)*y(n.VEGFR2_r) - p.koffVR*y(n.VEGFVEGFR2_r);

% Antibody and VEGF-Antibody complex
dydt(n.Ab_b)     = - p.konVA*y(n.Ab_b)*y(n.VEGF_b) + p.koffVA*y(n.VEGFAb_b) + (- p.k_bt*y(n.Ab_b)     - p.k_br*y(n.Ab_b)     + p.k_tb*y(n.Ab_t)*p.Vol_t/p.Vol_b     + p.k_rb*y(n.Ab_r)*p.Vol_r/p.Vol_b     )*p.AbEx - p.kcl_A*y(n.Ab_b);
dydt(n.Ab_t)     = - p.konVA*y(n.Ab_t)*y(n.VEGF_t) + p.koffVA*y(n.VEGFAb_t) + (+ p.k_bt*y(n.Ab_b)*p.Vol_b/p.Vol_t                            - p.k_tb*y(n.Ab_t)                            )*p.AbEx;
dydt(n.Ab_r)     = - p.konVA*y(n.Ab_r)*y(n.VEGF_r) + p.koffVA*y(n.VEGFAb_r) + (                       + p.k_br*y(n.Ab_b)*p.Vol_b/p.Vol_r                            - p.k_rb*y(n.Ab_r)     )*p.AbEx;

dydt(n.VEGFAb_b) = + p.konVA*y(n.Ab_b)*y(n.VEGF_b) - p.koffVA*y(n.VEGFAb_b) + (- p.k_bt*y(n.VEGFAb_b) - p.k_br*y(n.VEGFAb_b) + p.k_tb*y(n.VEGFAb_t)*p.Vol_t/p.Vol_b + p.k_rb*y(n.VEGFAb_r)*p.Vol_r/p.Vol_b )*p.AbEx - p.kcl_VA*y(n.VEGFAb_b);
dydt(n.VEGFAb_t) = + p.konVA*y(n.Ab_t)*y(n.VEGF_t) - p.koffVA*y(n.VEGFAb_t) + (+ p.k_bt*y(n.VEGFAb_b)*p.Vol_b/p.Vol_t                        - p.k_tb*y(n.VEGFAb_t)                        )*p.AbEx;
dydt(n.VEGFAb_r) = + p.konVA*y(n.Ab_r)*y(n.VEGF_r) - p.koffVA*y(n.VEGFAb_r) + (                       + p.k_br*y(n.VEGFAb_b)*p.Vol_b/p.Vol_r                        - p.k_rb*y(n.VEGFAb_r) )*p.AbEx;

% VEGF ligand
dydt(n.VEGF_b)   = - p.konVA*y(n.Ab_b)*y(n.VEGF_b) + p.koffVA*y(n.VEGFAb_b) - p.k_bt*y(n.VEGF_b)   - p.k_br*y(n.VEGF_b)   + p.k_tb*y(n.VEGF_t)*p.Vol_t/p.Vol_b   + p.k_rb*y(n.VEGF_r)*p.Vol_r/p.Vol_b   - p.kcl_V*y(n.VEGF_b) + p.VEGFprod_b - p.konVR*y(n.VEGF_b)*y(n.VEGFR2_b) + p.koffVR*y(n.VEGFVEGFR2_b);
dydt(n.VEGF_t)   = - p.konVA*y(n.Ab_t)*y(n.VEGF_t) + p.koffVA*y(n.VEGFAb_t) + p.k_bt*y(n.VEGF_b)*p.Vol_b/p.Vol_t                          - p.k_tb*y(n.VEGF_t)                                                + p.VEGFprod_t                 - p.konVR*y(n.VEGF_t)*y(n.VEGFR2_t) + p.koffVR*y(n.VEGFVEGFR2_t); 
dydt(n.VEGF_r)   = - p.konVA*y(n.Ab_r)*y(n.VEGF_r) + p.koffVA*y(n.VEGFAb_r)                        + p.k_br*y(n.VEGF_b)*p.Vol_b/p.Vol_r                          - p.k_rb*y(n.VEGF_r)                         + p.VEGFprod_r                 - p.konVR*y(n.VEGF_r)*y(n.VEGFR2_r) + p.koffVR*y(n.VEGFVEGFR2_r);

