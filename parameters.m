function parameters()
    % Define constants
    global pi MT_numb_M MT_numb_D Vg Vs Vs_c Pr_catastrophe Pr_rescue ...
           numberContactWindows numRegions regionAngles regionForceMultipliers ...
           regionProbabilities contactWindowAngles BasePosVecX BasePosVecY ...
           R1_max R2_max Prad F_MT Tau Duration Eta Mu D Fratio kM kD ...
           contact_length translation motherSpringOn daughterSpringOn;

    % Mathematical constants
    pi = 3.141592653589793;

    % MT parameters
    MT_numb_M = 30;  % Number of MTs (Mother)
    MT_numb_D = 30;  % Number of MTs (Daughter)
    Vg = 0.05;       % Growth velocity
    Vs = 0.005;      % Shrinkage velocity
    Vs_c = 0.002;    % Catastrophe shrinkage velocity
    Pr_catastrophe = 0.1;  % Probability of catastrophe
    Pr_rescue = 0.05;      % Probability of rescue

    % Contact parameters
    numberContactWindows = 5;
    numRegions = 6;
    regionAngles = [0.0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3, 2*pi];
    regionForceMultipliers = [1.0, 0.8, 1.2, 1.0, 0.8, 1.2];
    regionProbabilities = [0.3, 0.6, 0.9, 0.3, 0.6, 0.9];
    contactWindowAngles = [0.0, pi/2, pi, 3*pi/2, 2*pi, 5*pi/2];

    % Base position vectors
    BasePosVecX = 1.0;
    BasePosVecY = 1.0;

    % Geometric parameters
    R1_max = 2.0;  % Max radius along x
    R2_max = 2.0;  % Max radius along y
    Prad = 0.1;    % Pronucleus radius

    % Force and simulation parameters
    F_MT = 0.02;     % Force magnitude
    Tau = 0.01;      % Time step
    Duration = 100;  % Simulation duration
    Eta = 1.0;       % Drag coefficient
    Mu = 1.0;        % Viscosity coefficient
    D = 0.1;         % Diffusion coefficient
    Fratio = 1.0;    % Force ratio

    % Spring constants
    kM = 0.1;  % Spring constant for mother
    kD = 0.1;  % Spring constant for daughter
    contact_length = 1.0;

    % Simulation flags
    translation = true;
    motherSpringOn = true;
    daughterSpringOn = true;
    return 
end
