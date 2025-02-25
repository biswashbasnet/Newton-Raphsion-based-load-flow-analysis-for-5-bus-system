%% Biswash_Basnet Load Flow analysis  using Newton_Raphsion Method
clc;
clear all;
warning('off', 'all');

%% Predefined Line Data
BMva = 100;  % Base MVA
% Line data: {LineNo, StartBus, EndBus, R, X, Charging MVar}
predefined_data = { ...
    1, 'Birch', 'Elm',   0.042, 0.168, 4.1;
    2, 'Birch', 'Pine',  0.031, 0.126, 3.1;
    3, 'Elm',   'Maple', 0.031, 0.126, 3.1;
    4, 'Maple', 'Oak',   0.084, 0.336, 8.2;
    5, 'Maple', 'Pine',  0.053, 0.210, 5.1;
    6, 'Oak',   'Pine',  0.063, 0.252, 6.1; 
    };

% Convert predefined line data to table
lineData = cell2table(predefined_data, ...
    'VariableNames', {'LineNo', 'StartBus', 'EndBus', 'R', 'X', 'M'});

%% Process Line Data and Create Bus Mapping
% Extract start and end bus names from lineData
start_buses = lineData.StartBus;
end_buses   = lineData.EndBus;
unique_buses = unique([start_buses; end_buses]);% Get unique bus names in order
Nbus = length(unique_buses);
bus_map = containers.Map(unique_buses, 1:Nbus);% Create mapping: BusName --> Bus Number (1,2,...)
bus_reverse_map = containers.Map(1:Nbus, unique_buses);% For display, create a reverse map: Bus Number --> BusName

% Convert lineData to numeric matrix using the bus mapping
Nline = height(lineData);
linedata = zeros(Nline, 6);
for i = 1:Nline
    linedata(i,:) = [ ...
        lineData.LineNo(i), ...
        bus_map(char(lineData.StartBus{i})), ...
        bus_map(char(lineData.EndBus{i})), ...
        lineData.R(i), lineData.X(i), lineData.M(i) ];
end

%% Display Line Data Summary
disp(['Total Buses (from line data): ', num2str(Nbus)]);
disp(['Total Lines: ', num2str(Nline)]);
disp('---------------------------------------------------------------------------------------------------');
disp(' LineNo | Start Bus  | End Bus   | Resistance (p.u.) | Reactance (p.u.) | Charging MVar');
disp('---------------------------------------------------------------------------------------------------');
for row = 1:Nline
    sBus = linedata(row,2);
    eBus = linedata(row,3);
    fprintf(' %-6d | %-4d %-8s | %-4d %-8s | %-17.6f | %-17.6f | %-12.6f \n', ...
        linedata(row,1), sBus, bus_reverse_map(sBus), eBus, bus_reverse_map(eBus), ...
        linedata(row,4), linedata(row,5), linedata(row,6));
end
disp('---------------------------------------------------------------------------------------------------');

%% Ybus Matrix Calculation (Including Charging Mvar)
Ybus = zeros(Nbus, Nbus); 
for i = 1:Nline
    p = linedata(i, 2);
    q = linedata(i, 3);
    Z_line = linedata(i, 4) + 1j * linedata(i, 5);
    if Z_line ~= 0  % Avoid division by zero
        yline = 1 / Z_line ;
        Ybus(p, p) = Ybus(p, p) + yline;
        Ybus(q, q) = Ybus(q, q) + yline;
        Ybus(p, q) = Ybus(p, q) - yline;
        Ybus(q, p) = Ybus(q, p) - yline;
    end
end

%Charging MVar (half at each bus)
Y_shunt =  zeros(Nbus, Nbus);
for k = 1:Nline
    p = linedata(k, 2);
    q = linedata(k, 3);
   y_shunt = 1j * linedata(k, 6)/BMva;
 
    Y_shunt(p, p)= Y_shunt(p, p)  + y_shunt/2 ;
   Y_shunt(q, q) = Y_shunt(q, q) + y_shunt/2 ;
end

%add charging MVar, final bus
Ybus = Ybus + Y_shunt;
B = zeros(Nbus,Nbus);

for k = 1:Nline
    p = linedata(k, 2);
    q = linedata(k, 3);
    b = 1j * linedata(k, 6)/BMva;
    B(p,q)= b/2;
    B(q,p)= b/2;
end

% Compute magnitude and angle of Ybus for later use
Ymag = abs(Ybus);
theta = angle(Ybus);
disp(' Ybus matrix incorporating charging MVar:');
disp(Ybus);

%% Bus Data Definition with Bus Names
% Bus data: {BusName, BusType, Pg, Qg, Pd, Qd, |V|, delta (deg), Qmin, Qmax}
busdata_cell = { ... 
    'Elm',    3, 0,    0,  1.15,  0.6,  1.00, 0,  0,  0;
    'Oak',    3, 0,    0,  0.7,   0.3,  1.00, 0,  0,  0;
    'Maple',  2, 1.8,  0,  0.7,   0.4,  1.02, 0, -2,  2;
    'Birch',  1, 0,    0,  0.65,  0.30, 1.04, 0,  0,  0;
    'Pine',   3, 0,    0,  0.85,  0.4,  1.00, 0,  0,  0;
    };

% Bus Type: Slack = 1, PV = 2, PQ = 3
busdata = cell2table(busdata_cell, 'VariableNames', ...
    {'BusName','BusType','Pg','Qg','Pd','Qd','Vmag','delta_deg','Qmin','Qmax'});

%% Reorder Bus Data According to Ybus Mapping
ordered_busdata = table();
for i = 1:Nbus
    currentBusName = bus_reverse_map(i);
    idx = find(strcmp(busdata.BusName, currentBusName));
    if isempty(idx)
        error('Bus name %s not found in busdata.', currentBusName);
    end
    ordered_busdata = [ordered_busdata; busdata(idx,:)];
end
busdata = ordered_busdata;
disp('Ordered Bus Data:');
disp(busdata);

%% Bus Data Initialization for Load Flow
type = busdata.BusType;
Pg   = busdata.Pg;    % Generation (p.u.)
Qg   = busdata.Qg;
Pd   = busdata.Pd;    % Load (p.u.)
Qd   = busdata.Qd;
Qmin = busdata.Qmin;
Qmax = busdata.Qmax;
Vmag = busdata.Vmag;  % Voltage magnitudes (p.u.)
delta = deg2rad(busdata.delta_deg);  % Voltage angles (radians)
V = Vmag .* (cos(delta) + 1i*sin(delta));

% Scheduled net injections (Generation minus Load)
Psch = Pg - Pd;
Qsch = Qg - Qd;

%% Newton-Raphson Load Flow Iteration Parameters
accuracy = 1;
iter = 1;
max_iter = 10;  % Maximum number of iterations

while accuracy >= 1e-2 && iter <= max_iter
    fprintf('\n========== Iteration %d ==========\n', iter);
    
    %% Calculate Bus Injections (Calculated)
    P_cal = zeros(Nbus, 1);
    Q_cal = zeros(Nbus, 1);
    for i = 1:Nbus
        for n = 1:Nbus
            P_cal(i) = P_cal(i) + Vmag(i)*Vmag(n)*Ymag(i,n)* ...
                       cos(theta(i,n) + delta(n) - delta(i));
            Q_cal(i) = Q_cal(i) - Vmag(i)*Vmag(n)*Ymag(i,n)* ...
                       sin(theta(i,n) + delta(n) - delta(i));
        end
    end
    
    %% Q Limit Checking for PV Buses
    for i = 1:Nbus
        if type(i) == 2  % PV bus
            if Q_cal(i) > Qmax(i)
                Q_cal(i) = Qmax(i);
                type(i) = 3;  % Convert to PQ temporarily if exceeded
            elseif Q_cal(i) < Qmin(i)
                Q_cal(i) = Qmin(i);
                type(i) = 3;
            else
                type(i) = 2;  % Remains PV
            end
        end
    end
    
    %% Identify Non-Slack and PQ Buses
    non_slack = find(type ~= 1);   % PV and PQ buses
    pq = find(type == 3);          % PQ buses only
    n_non_slack = length(non_slack);
    npq = length(pq);
    
    %% Calculate Mismatch Vectors (ΔP and ΔQ)
    DP = Psch(non_slack) - P_cal(non_slack);
    DQ = Qsch(pq) - Q_cal(pq);
    mismatch = [DP; DQ]; % disp(mismatch);
    
    %% Jacobian Matrix Calculation
    % J1: dP/dδ for non-slack buses
    J1 = zeros(n_non_slack, n_non_slack);
    for i = 1:n_non_slack
        ibus = non_slack(i);
        for j = 1:n_non_slack
            jbus = non_slack(j);
            if ibus == jbus
                J1(i,j) = -Q_cal(ibus) - Vmag(ibus)^2 * imag(Ybus(ibus,ibus));
            else
                J1(i,j) = -Vmag(ibus)*Vmag(jbus)*Ymag(ibus,jbus)* ...
                          sin(theta(ibus,jbus) + delta(jbus) - delta(ibus));
            end
        end
    end
    
    % J2: dP/dV for non-slack buses with respect to PQ buses
    J2 = zeros(n_non_slack, npq);
    for i = 1:n_non_slack
        ibus = non_slack(i);
        for j = 1:npq
            jbus = pq(j);
            if ibus == jbus
                J2(i,j) = (P_cal(ibus)/Vmag(ibus)) + Vmag(ibus)*real(Ybus(ibus,ibus));
            else
                J2(i,j) = Vmag(ibus)*Ymag(ibus,jbus)* ...
                          cos(theta(ibus,jbus) + delta(jbus) - delta(ibus));
            end
        end
    end
    
    % J3: dQ/dδ for PQ buses with respect to nos slack bus
    J3 = zeros(npq, n_non_slack);
    for i = 1:npq
        ibus = pq(i);
        for j = 1:n_non_slack
            jbus = non_slack(j);
            if ibus == jbus
                J3(i,j) = P_cal(ibus) - Vmag(ibus)^2 * real(Ybus(ibus,ibus));
            else
                J3(i,j) = -Vmag(ibus)*Vmag(jbus)*Ymag(ibus,jbus)* ...
                          cos(theta(ibus,jbus) + delta(jbus) - delta(ibus));
            end
        end
    end
    
    % J4: dQ/dV for PQ buses
    J4 = zeros(npq, npq);
    for i = 1:npq
        ibus = pq(i);
        for j = 1:npq
            jbus = pq(j);
            if ibus == jbus
                J4(i,j) = (Q_cal(ibus)/Vmag(ibus)) - Vmag(ibus)*imag(Ybus(ibus,ibus));
            else
                J4(i,j) = -Vmag(ibus)*Ymag(ibus,jbus)* ...
                          sin(theta(ibus,jbus) + delta(jbus) - delta(ibus));
            end
        end
    end
    
    % Form the full Jacobian matrix
    J = [J1, J2; J3, J4];
    disp('Jacobian matrix:');
    disp(J);
    
    %% Solve for State Variable Updates (Δδ and ΔV)
    DX = J \ mismatch;
    dDelta = DX(1:n_non_slack);
    dV = DX(n_non_slack+1:end);
    
    %% Update State Variables
    delta(non_slack) = delta(non_slack) + dDelta;
    Vmag(pq) = Vmag(pq) + dV;
    V = Vmag .* (cos(delta) + 1i*sin(delta));
    
    %% Check Convergence
    accuracy = max(abs(mismatch));
    fprintf('Max mismatch = %g\n', accuracy);
    iter = iter + 1;
end

%% Display Final Bus Voltages
Del = rad2deg(delta);  % Convert radians to degrees
disp('Final Load Flow Solution (Newton-Raphson):');
disp('Bus   | V (pu)  | Angle (deg)');
for i = 1:Nbus
    fprintf('%2d    %7.4f   %7.4f\n', i, Vmag(i), Del(i));
end

%% Bus Power Injections, Generation, and Loads (p.u.)
% Compute bus injections: S = V .* conj(I)
Ibus = Ybus * V;
Si = V .* conj(Ibus) * BMva;  % In MVA
Pi = real(Si);
Qi = imag(Si);

% Scheduled load (in MW and MVar)
Pload_disp = Pd * BMva;
Qload_disp = Qd * BMva;

% Actual generation (Injection + Load; Pgen = P_injection + P_load)
Pgen_actual = Pi + Pload_disp;
Qgen_actual = Qi + Qload_disp;

disp('---------------------------------------------------------------------------------------------------------------------------');
disp('                              Newton-Raphson based Load Flow Analysis');
disp('--------------------------------------------------------------------------------------------------------------------------');
disp('| Bus |  V (pu)  | Angle (deg) | Injection (MW) | Injection (MVar) | Gen (MW) | Gen (MVar) | Load (MW) | Load (MVar) |');
for i = 1:Nbus
    fprintf('%3d    %8.4f    %8.4f      %8.3f         %8.3f         %8.3f   %8.3f    %8.3f     %8.3f\n', ...
        i, Vmag(i), Del(i), Pi(i), Qi(i), Pgen_actual(i), Qgen_actual(i), Pload_disp(i), Qload_disp(i));
end

%% Line Flow Calculations
% Calculate branch flows using: I_pq = -(V(p)-V(q))*Ybus(p,q)
start_no = linedata(:,2);
end_no = linedata(:,3);
Sij = zeros(Nbus, Nbus);  % Branch flow matrix in MVA

for m = 1:Nline
    p = start_no(m);
    q = end_no(m);
    I(p,q) = -(V(p) - V(q)) * Ybus(p,q);% Since Ybus(p,q) = -y_line for a branch between p and q, 
    I(q,p) = -I(p,q);
    
end
for m=1: Nbus
    for n= 1: Nbus
        if m~=n 

          I(m,n) =   I(m,n)+ V(m)* B(m,n);
          Powerflow(m,n) = V(m) * conj(I(m,n)) * BMva;
        end
    end

end

% Compute line losses for each branch: Loss = Powerflow(p,q) + Powerflow(q,p)
Lij = zeros(Nline,1);
for m = 1:Nline
    p = start_no(m);
    q = end_no(m);
    Lij(m) = Powerflow(p,q) + Powerflow(q,p);
end
Lpij = real(Lij);
Lqij = imag(Lij);

disp('----------------------------------------------------------------------------------------------------------------------');
disp('                              Line Flow and Losses Calculation');
disp('-----------------------------------------------------------------------------------------------------------------------');
disp('| From Bus | To Bus | Flow (MW) | Flow (MVar) | Reverse Flow (MW) | Reverse Flow (MVar) | Loss (MW) | Loss (MVar) |');
for m = 1:Nline
    p = start_no(m);
    q = end_no(m);
    fprintf('%9d  %7d   %10.3f    %10.3f       %10.3f         %10.3f       %10.3f   %10.3f\n', ...
        p, q, real(Powerflow(p,q)), imag(Powerflow(p,q)), real(Powerflow(q,p)), imag(Powerflow(q,p)), Lpij(m), Lqij(m));
end
disp('---------------------------------------------------------------------------------------------------------------------');
fprintf('   Total Loss:    %10.3f MW    %10.3f MVar\n', sum(Lpij), sum(Lqij));
disp('########################################################################################################################');
