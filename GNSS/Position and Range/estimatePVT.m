function [P,V,T,DOP,sigma,N_sats,fit] = estimatePVT(x_0,b_0,x,T_0,nav,obs,t,stop_condition,opts)
arguments
    x_0 (1,3) double
    b_0 (1,1) double
    x (1,3) double
    T_0 (1,1) double
    nav (1,1) struct
    obs (1,1) struct
    t (:,2) double
    stop_condition (1,1) string {mustBeMember(stop_condition,{'convergence','fixed'})} = 'convergence'
    opts.AbsTol (1,1) double = 1e-5
    opts.NumIter (1,1) double = 5
end

N_EPOCHS = size(t,1);

el_cutoff = 10;

P = nan(N_EPOCHS,3);
V = nan(N_EPOCHS,3);
T = nan(N_EPOCHS,1);
DOP = struct('P',zeros(N_EPOCHS,1),'T',zeros(N_EPOCHS,1),'G',zeros(N_EPOCHS,1),'E',zeros(N_EPOCHS,1),'N',zeros(N_EPOCHS,1),'V',zeros(N_EPOCHS,1),'H',zeros(N_EPOCHS,1));
sigma = zeros(N_EPOCHS,1);
N_sats = zeros(N_EPOCHS,1);
fit = struct('pre',[],'post',[],'el',[],'t',[]);

R_L = rot_enu(x);
R_L_aug = [R_L, zeros(3,1); zeros(1,3), 1];

for i = 1:N_EPOCHS
    t_ind = find((obs.TOW == t(i,1) & obs.WN == t(i,2)) & (obs.PRN < 100 & (obs.C1 ~= 0 & obs.P2 ~= 0)));

    PRN = obs.PRN(t_ind);
    C1 = obs.C1(t_ind);
    P2 = obs.P2(t_ind);
    
    R = zeros(length(t_ind),1);
    x_s_r = zeros(length(t_ind),3);
    x_s_t = zeros(length(t_ind),3);
    b_sv = zeros(length(t_ind),1);
    rel_sv = zeros(length(t_ind),1);
    
    x_hat = x_0;
    b_hat = b_0;

    n = 0;
    cond = true;
    while cond
        for tt = 1:length(t_ind)
            [x_s_r(tt,:),~,b_sv(tt),rel_sv(tt)] = eph2pos(nav,t(i,:),PRN(tt));
            [R(tt),x_s_t(tt,:)] = computeExpectedRange(x_hat,x_s_t(tt,:),nav,t(i,:),PRN(tt));
        end
    
        [~,el,~] = computeAzElRange(x_hat,x_s_t);

        ind_rem = el < el_cutoff;

        PRN(ind_rem) = [];
        C1(ind_rem) = [];
        P2(ind_rem) = [];
        R(ind_rem) = [];
        x_s_t(ind_rem,:) = [];
        b_sv(ind_rem) = [];
        rel_sv(ind_rem) = [];
        t_ind(ind_rem) = [];
        el(ind_rem) = [];

        [rho_IF,~,~,~] = ionoModel(C1,f_C('gps','L1'),P2,f_C('gps','L2'));
        T_err = tropoModel(T_0,el);

        % need negative sign in front of position diff, shouldn't include clock bias here?
        A = [(x_s_t - x_hat)./(R - b_hat),ones(length(t_ind),1)];
        delta_rho_pre = rho_IF - (R - b_sv - rel_sv + T_err) - b_hat;
        delta_xb_hat = (A.'*A)\(A.'*delta_rho_pre);

        delta_x_hat = delta_xb_hat(1:3).';
        delta_b_hat = delta_xb_hat(4);
    
        x_hat = x_hat + delta_x_hat;
        b_hat = b_hat + delta_b_hat;
    
        switch stop_condition
            case 'convergence'
                res = norm(delta_xb_hat);
                cond = res > opts.AbsTol;
            case 'fixed'
                n = n + 1;
                cond = n < opts.NumIter;
        end
    end

    delta_rho_post = delta_rho_pre - A*delta_xb_hat;
    sigma_rho = std(delta_rho_post);

    Q = inv(A.'*A);

    PDOP = sqrt(Q(1,1) + Q(2,2) + Q(3,3));
    TDOP = sqrt(Q(4,4));
    GDOP = sqrt(sum(diag(Q)));

    A_enu = A*R_L_aug.';
    Q_enu = inv(A_enu.'*A_enu);

    EDOP = sqrt(Q_enu(1,1));
    NDOP = sqrt(Q_enu(2,2));
    VDOP = sqrt(Q_enu(3,3));
    HDOP = sqrt(Q_enu(1,1) + Q_enu(2,2));

    P(i,:) = x_hat;
    V(i,:) = [0,0,0];
    T(i) = b_hat;
    DOP.P(i) = PDOP;
    DOP.T(i) = TDOP;
    DOP.G(i) = GDOP;
    DOP.E(i) = EDOP;
    DOP.N(i) = NDOP;
    DOP.V(i) = VDOP;
    DOP.H(i) = HDOP;
    sigma(i) = sigma_rho;
    N_sats(i) = length(t_ind);

    fit.pre = [fit.pre; delta_rho_pre];
    fit.post = [fit.post; delta_rho_post];
    fit.el = [fit.el; el];
    fit.t = [fit.t; t(i,1)*ones(length(t_ind),1)];
end

end
