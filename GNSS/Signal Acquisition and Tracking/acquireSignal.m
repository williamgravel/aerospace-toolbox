function [f_D,tau,S] = acquireSignal(s,PRN,T,f_s,f_IF,opts)
arguments
    s (1,:) double
    PRN (1,1) double
    T (1,1) double
    f_s (1,1) double
    f_IF (1,1) double
    opts.FreqSearchRange (1,2) double = [-10e3,10e3]
    opts.FreqSearchStep (1,1) double = 1/(2*T);
    opts.FreqType (1,1) string {mustBeMember(opts.FreqType,["doppler","total"])} = "doppler"
    opts.CodeSearchRange (1,2) double = [0,T*1.023e6]
    opts.CodeSearchStep (1,1) double = 1/2
    opts.CodeType (1,1) string {mustBeMember(opts.CodeType,["chip","sample","time"])} = "chip"
    opts.SearchType (1,1) string {mustBeMember(opts.SearchType,["serial","parallel"])} = "parallel"
    opts.SignalType (1,1) string {mustBeMember(opts.SignalType,["real","complex","conjugate"])} = "conjugate"
    opts.CorrelatorType (1,1) string {mustBeMember(opts.CorrelatorType,["raw","magnitude","normalized"])} = "normalized"
    opts.PlotCorrelator (1,1) logical = false
    opts.PlotCorrelatorType (1,1) string {mustBeMember(opts.PlotCorrelatorType,["bar","mesh","surf"])} = "mesh"
    opts.PlotPeakSlices (1,1) logical = false
end

% Manipulate input signal if requested
switch opts.SignalType
    case "real"
        s = real(s);
    case "complex"
        % do nothing
    case "conjugate"
        s = conj(s);
end

% Define chipping rate and time step sizes
f_c = 1.023e6; % chipping rate [Hz]
Delta_t_c = 1/f_c; % chip duration [s]
Delta_t_s = 1/f_s; % sample duration [s]

% Define number of samples and time vectors
N = T*f_s; % number of samples
t = 0:Delta_t_s:(T - Delta_t_s); % sampling time vector [s]
t_c = 0:Delta_t_c:(T - Delta_t_c); % chipping time vector [s]

% Generate and interpolate PRN C/A code sequence
x = generateCACode(PRN,length(t_c),'CodeType','signed');
x = interp1(t_c,x,t,'previous','extrap');

% tic
if strcmp(opts.SearchType,"serial")
    % Construct Doppler shift and chip delay mesh grids
    f_vec = opts.FreqSearchRange(1):opts.FreqSearchStep:opts.FreqSearchRange(2);
    d_vec = opts.CodeSearchRange(1):opts.CodeSearchStep:opts.CodeSearchRange(2);
    [D,F] = meshgrid(floor(d_vec*Delta_t_c/Delta_t_s),f_vec + f_IF);
    
    % Initialize autocorrelation matrix
    S = zeros(length(f_vec),length(d_vec));
    
    % Compute complex angular frequency grid
    OMEGA = -1j*2*pi*F; % angular frequency grid [rad/s]
    
    % Compute autocorrelation values for each grid point
    for i = 1:N
        S = S + s(i + D).*x(i).*exp(OMEGA*t(i));
    end
elseif strcmp(opts.SearchType,"parallel")
    % Construct Doppler shift and chip delay mesh grids
    f_vec = (0:N-1)*(f_s/N) - f_IF;
    d_vec = opts.CodeSearchRange(1):opts.CodeSearchStep:opts.CodeSearchRange(2);
    [D,~] = meshgrid(floor(d_vec*Delta_t_c/Delta_t_s),x);

    % Compute FFT of delayed signals
    S = fft(s((1:N)' + D).*x');

    % Truncate frequency bins outside of search range
    S(f_vec < opts.FreqSearchRange(1) | f_vec > opts.FreqSearchRange(2),:) = [];
    f_vec(f_vec < opts.FreqSearchRange(1) | f_vec > opts.FreqSearchRange(2)) = [];
end
% toc

% Determine corresponding Doppler shift and chip delay
[~,I] = max(abs(S),[],'all');
[I_row,I_col] = ind2sub(size(S),I);
f_D = f_vec(I_row);
tau = d_vec(I_col);

switch opts.CodeType
    case "chip"
        % do nothing
    case "sample"
        tau = floor(tau*Delta_t_c/Delta_t_s);
    case "time"
        tau = tau*Delta_t_c;
end

% Manipulate autocorrelation values
switch opts.CorrelatorType
    case "raw"
        % do nothing
    case "magnitude"
        S = abs(S);
    case "normalized"
        S = abs(S)/abs(S(I));
end

% Plot results
if opts.PlotCorrelator
    figure('Units','normalized','Position',[1/4 1/4 1/2 1/2])
    
    switch opts.PlotCorrelatorType
        case "bar"
            b = bar3(f_vec,abs(S));
            set(b,{'CData'}',get(b,'ZData'));
            set(b,'FaceColor','interp');
            colorbar
        case "mesh"
            switch opts.CodeType
                case "chip"
                    mesh(d_vec,f_vec,abs(S));
                    xlim(opts.CodeSearchRange)
                case "sample"
                    mesh(floor(d_vec*Delta_t_c/Delta_t_s),f_vec,abs(S));
                    xlim(floor(opts.CodeSearchRange*Delta_t_c/Delta_t_s))
                case "time"
        
            end
            hidden off
            colorbar
        case "surf"
            surf(d_vec,f_vec,abs(S));
            colorbar
    end
    xlabel(sprintf('code delay [%s]',opts.CodeType))
    ylabel('Doppler shift [Hz]')
    zlabel('autocorrelation')
    ylim(opts.FreqSearchRange)
    pbaspect([1,1,0.5])
    view(60,15)
    title(sprintf('C/A CODE PRN %2d %.0f MS',PRN,T*1000),sprintf('PEAK @ %.1f %ss, %.0f Hz',tau,opts.CodeType,f_D))
end

if opts.PlotPeakSlices
    figure('Units','normalized','Position',[1/4 1/4 1/2 1/2])
    sgtitle(sprintf('C/A CODE PRN %2d %.0f MS',PRN,T*1000))

    subplot(2,1,1)
    switch opts.CodeType
        case "chip"
            plot(d_vec,S(I_row,:),'b','LineWidth',1)
            xlim(opts.CodeSearchRange)
        case "sample"
            plot(floor(d_vec*Delta_t_c/Delta_t_s),S(I_row,:),'b','LineWidth',1)
            xlim(floor(opts.CodeSearchRange*Delta_t_c/Delta_t_s))
        case "time"

    end
    grid on
    ylim([0,1])
    xlabel(sprintf('code delay [%s]',opts.CodeType))
    ylabel('autocorrelation')
    title('CONSTANT DOPPLER SHIFT SLICE',sprintf('PEAK @ %.0f Hz',f_D))

    subplot(2,1,2)
    plot(f_vec,S(:,I_col),'r','LineWidth',1)
    grid on
    xlim(opts.FreqSearchRange)
    ylim([0,1])
    xlabel('Doppler shift [Hz]')
    ylabel('autocorrelation')
    title('CONSTANT CODE DELAY SLICE',sprintf('PEAK @ %.1f %ss',tau,opts.CodeType))
end

end
