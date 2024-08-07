function main()

close all;
clear all;
    % Main function to run the simulation
    global MT_Pos_M MT_Pos_D force_M force_D proNucPos psi force torque ...
           torque_M torque_D basePosM basePosD;

    global pi MT_numb_M MT_numb_D Vg Vs Vs_c Pr_catastrophe Pr_rescue ...
           numberContactWindows numRegions regionAngles regionForceMultipliers ...
           regionProbabilities contactWindowAngles BasePosVecX BasePosVecY ...
           R1_max R2_max Prad F_MT Tau Duration Eta Mu D Fratio kM kD ...
           contact_length translation motherSpringOn daughterSpringOn;
    
    global MT_Pos_M MT_Pos_D MT_Growing_M MT_Growing_D MT_GrowthVel_M ...
           MT_GrowthVel_D MT_Contact_M MT_Contact_D basePosM basePosD ...
           force_M force_D proNucPos psi force torque torque_M torque_D contacts ...
           regionPos generator stdNormalDist;

    initialization(100,100);  % Initialize parameters and variables

    %% I tried with 10000,10000 and led to computer hanging
    
    %% Simulation settings
    numRuns = 1;      % Number of simulation runs
    writeAllData = false; % Flag to write all data
    writeTempData = false; % Flag to write temporary data
    ONLY_COMMA = true; % Output formatting flag
    fileName = 'output.txt'; % Output file name

    % Open file for writing data
    fileID = fopen(fileName, 'w');

    for run = 1:numRuns
        % Run the model simulation
        runModel(writeAllData, writeTempData, fileID, ONLY_COMMA);

        % Reset to base positions after each run
        setToBasePos();
    end

    fclose(fileID); % Close the output file
end

function runModel(writeAllData, writeTempData, fileID, ONLY_COMMA)
    % Run the model simulation
    global MT_Pos_M MT_Pos_D force_M force_D proNucPos psi force torque ...
           torque_M torque_D basePosM basePosD Tau Duration Prad;

    nextWrite = 0;
    writeInterval = 0.5;

    if writeAllData
        % Write full data if specified
        writeData(0, fileID, ONLY_COMMA);
    elseif writeTempData
        % Write partial data and schedule next write
        writePartialData(0, fileID, ONLY_COMMA);
        nextWrite = nextWrite + writeInterval;
    end

    % Simulation loop
    for t = 0:Tau:Duration
        % Update the pronuclear position
        updatePNPos();

        % Handle MT growth and shrinkage
        handleMTGrowth();

        % Write data at specified intervals
        if writeAllData
            writeData(t + Tau, fileID, ONLY_COMMA);
        elseif writeTempData && (t > nextWrite)
            writePartialData(t + Tau, fileID, ONLY_COMMA);
            nextWrite = nextWrite + writeInterval;
        end
    end

    % Write final data if required
    if ~(writeAllData || writeTempData)
        writeFinalData(fileID, ONLY_COMMA);
    end
end

function writeData(t, fileID, ONLY_COMMA)
    % Write all relevant data to the file
    global MT_Pos_M MT_Pos_D force_M force_D proNucPos psi force torque ...
           torque_M torque_D basePosM basePosD;

    fprintf(fileID, '%f', t);
    fprintf(fileID, formatOutput(proNucPos, ONLY_COMMA));
    fprintf(fileID, formatOutput(psi, ONLY_COMMA));
    fprintf(fileID, formatOutput(MT_Pos_M, ONLY_COMMA));
    fprintf(fileID, formatOutput(MT_Pos_D, ONLY_COMMA));
    fprintf(fileID, formatOutput(force_M, ONLY_COMMA));
    fprintf(fileID, formatOutput(force_D, ONLY_COMMA));
    fprintf(fileID, formatOutput(force, ONLY_COMMA));
    fprintf(fileID, formatOutput(torque_M, ONLY_COMMA));
    fprintf(fileID, formatOutput(torque_D, ONLY_COMMA));
    fprintf(fileID, formatOutput(torque, ONLY_COMMA));
    fprintf(fileID, formatOutput(basePosM, ONLY_COMMA));
    fprintf(fileID, formatOutput(basePosD, ONLY_COMMA));
    fprintf(fileID, '\n');
end

function writeFinalData(fileID, ONLY_COMMA)
    % Write final position data
    global proNucPos psi;
    fprintf(fileID, formatOutput(proNucPos, ONLY_COMMA));
    fprintf(fileID, formatOutput(psi, ONLY_COMMA));
    fprintf(fileID, '\n');
end

function writePartialData(t, fileID, ONLY_COMMA)
    % Write some relevant positional data
    global proNucPos psi basePosM basePosD;
    fprintf(fileID, '%f', t);
    fprintf(fileID, formatOutput(proNucPos, ONLY_COMMA));
    fprintf(fileID, formatOutput(psi, ONLY_COMMA));
    fprintf(fileID, formatOutput(basePosM, ONLY_COMMA));
    fprintf(fileID, formatOutput(basePosD, ONLY_COMMA));
    fprintf(fileID, '\n');
end

function mtForceCalcM(forceMag)
    % Calculate forces on MTOC M due to MTs
    global MT_Pos_M basePosM force_M MT_Contact_M regionAngles regionForceMultipliers numRegions;

    force_M = [0.0, 0.0];
    for i = 1:length(MT_Contact_M)
        if MT_Contact_M(i) == 0
            continue;
        end

        % Calculate angle and determine force multiplier
        angle = atan2(MT_Pos_M(i, 2), MT_Pos_M(i, 1));
        if angle < 0
            angle = angle + 2 * pi;
        end

        multiplier = 1;
        for j = 2:numRegions
            if angle < regionAngles(j) && angle >= regionAngles(j-1)
                multiplier = regionForceMultipliers(j-1);
            end
        end

        % Calculate force
        mtUnitVec = (MT_Pos_M(i, :) - basePosM) / norm(MT_Pos_M(i, :) - basePosM);
        force_M = force_M + multiplier * mtUnitVec;
    end
    force_M = force_M * forceMag;
end

function mtForceCalcD(forceMag)
    % Calculate forces on MTOC D due to MTs
    global MT_Pos_D basePosD force_D MT_Contact_D regionAngles regionForceMultipliers numRegions;

    force_D = [0.0, 0.0];
    for i = 1:length(MT_Contact_D)
        if MT_Contact_D(i) == 0
            continue;
        end

        % Calculate angle and determine force multiplier
        angle = atan2(MT_Pos_D(i, 2), MT_Pos_D(i, 1));
        if angle < 0
            angle = angle + 2 * pi;
        end

        multiplier = 1;
        for j = 2:numRegions
            if angle < regionAngles(j) && angle >= regionAngles(j-1)
                multiplier = regionForceMultipliers(j-1);
            end
        end

        % Calculate force
        mtUnitVec = (MT_Pos_D(i, :) - basePosD) / norm(MT_Pos_D(i, :) - basePosD);
        force_D = force_D + multiplier * mtUnitVec;
    end
    force_D = force_D * forceMag;
end

function netMTForce(centrosome)
    % Wrapper to compute the net force on an MTOC
    global Fratio F_MT;

    switch centrosome
        case 'M'
            mtForceCalcM(Fratio * F_MT);
        case 'D'
            mtForceCalcD(F_MT);
    end
end

function updatePNPos()
    % Update the pronuclear position
    global proNucPos psi force torque force_M force_D torque_M torque_D ...
           basePosM basePosD motherSpringOn daughterSpringOn kM kD Eta Eta2 ...
           Tau D MT_Pos_M MT_Pos_D contacts regionAngles regionForceMultipliers ...
           Prad translation Mu;

    proNucRad = [Prad * cos(psi), Prad * sin(psi)];
    force = [0.0, 0.0];

    % Compute net forces on the MTOCs
    netMTForce('M');
    netMTForce('D');

    % Update positions for MTOC M
    if motherSpringOn
        springAnchorM = proNucPos + proNucRad;
        springForceM = kM * (springAnchorM - basePosM);
        force_M = force_M + springForceM;
        force = force - springForceM;
        % basePosM = basePosM + (1.0/Eta2) * Tau * force_M + sqrt(2 * D * Tau) * randn(1, 2);
        basePosM = basePosM + (1.0/Eta) * Tau * force_M + sqrt(2 * D * Tau) * randn(1, 2);

        torque_M = springForceM(1) * proNucRad(2) - springForceM(2) * proNucRad(1);
    else
        basePosM = proNucPos + proNucRad;
        force = force + force_M;
        torque_M = -force_M(1) * proNucRad(2) + force_M(2) * proNucRad(1);
    end

    % Update positions for MTOC D
    if daughterSpringOn
        springAnchorD = proNucPos - proNucRad;
        springForceD = kD * (springAnchorD - basePosD);
        force_D = force_D + springForceD;
        force = force - springForceD;
        % basePosD = basePosD + (1.0/Eta2) * Tau * force_D + sqrt(2 * D * Tau) * randn(1, 2);
        basePosD = basePosD + (1.0/Eta) * Tau * force_D + sqrt(2 * D * Tau) * randn(1, 2);

        torque_D = -springForceD(1) * proNucRad(2) + springForceD(2) * proNucRad(1);
    else
        basePosD = proNucPos - proNucRad;
        force = force + force_D;
        torque_D = force_D(1) * proNucRad(2) - force_D(2) * proNucRad(1);
    end
% % Here's the continuation and completion of the MATLAB translation for the `main.cpp` file. This will cover the entire functionality required for the simulation process:
% % 
% % main.m` (Continued matlab

    % Calculate displacements and torques
    if translation
        proNucPos = proNucPos + (1.0/Eta) * Tau * force;
    end

    torque = torque_M + torque_D;
    psi = psi + (1/Mu) * torque * Tau;
end

function advanceMT(vel, vec, mag)
global Tau, vel, mag, vec 
    % Advance an MT by growing or shrinking
    vec = vec * (1 + vel * Tau / mag);
end

function prob = probContact(ang)
    % Compute the probability of making contact at angle ang
    global contactWindowAngles regionProbabilities contacts numRegions;

    % Check contact windows
    for i = 2:length(contactWindowAngles)
        if ang < contactWindowAngles(i)
            if contacts(i-1)
                prob = 0;
                return;
            else
                break;
            end
        end
    end

    % Check regions
    for i = 2:numRegions
        if ang < regionAngles(i)
            prob = regionProbabilities(i-1);
            return;
        end
    end

    % Fail-safe return
    prob = 0;
end

function addContact(angle)
    % Add a contact at the specified angle
    global contacts contactWindowAngles;

    for i = 2:length(contactWindowAngles)
        if angle < contactWindowAngles(i)
            assert(~contacts(i-1));  % Assert no contact has been made
            contacts(i-1) = true;
            break;
        end
    end
end

function removeContact(angle)
    % Remove a contact at the specified angle
    global contacts contactWindowAngles;

    for i = 2:length(contactWindowAngles)
        if angle < contactWindowAngles(i)
            assert(contacts(i-1));  % Assert contact exists
            contacts(i-1) = false;
            break;
        end
    end
end

function mtContactTest(centrosome, i)
    % Test MT contact for the ith MT from the specified centrosome
    global MT_Pos_M MT_Pos_D MT_Contact_M MT_Contact_D MT_Growing_M ...
           MT_Growing_D MT_GrowthVel_M MT_GrowthVel_D contact_length Vs_c;

    switch centrosome
        case 'M'
            angleM = atan2(MT_Pos_M(i, 2), MT_Pos_M(i, 1));
            if angleM < 0
                angleM = angleM + 2 * pi;
            end

            if rand < probContact(angleM)
                MT_GrowthVel_M(i) = 0;
                MT_Contact_M(i) = contact_length;
                MT_Growing_M(i) = false;
                addContact(angleM);
            else
                MT_GrowthVel_M(i) = -Vs_c;
                MT_Growing_M(i) = false;
            end

        case 'D'
            angleD = atan2(MT_Pos_D(i, 2), MT_Pos_D(i, 1));
            if angleD < 0
                angleD = angleD + 2 * pi;
            end

            if rand < probContact(angleD)
                MT_GrowthVel_D(i) = 0;
                MT_Contact_D(i) = contact_length;
                MT_Growing_D(i) = false;
                addContact(angleD);
            else
                MT_GrowthVel_D(i) = -Vs_c;
                MT_Growing_D(i) = false;
            end
    end
end

function respawnMTB(vec, ang, envelope)
    % Base respawn function for an MT
    randR = rand;
    randT = rand;
    r = sqrt(randR);
    t = (envelope(2) - envelope(1)) * randT + ang + envelope(1) - pi/2.0;
    vec = [r * cos(t), r * sin(t)];
end

function respawnMT(centrosome, vec, i, envelope)
    % Wrapper function for respawnMTB for a specific MT
    global MT_Growing_M MT_Growing_D MT_GrowthVel_M MT_GrowthVel_D Vg;

    switch centrosome
        case 'M'
            respawnMTB(vec, psi, envelope);
            MT_Growing_M(i) = true;
            MT_GrowthVel_M(i) = Vg;
            MT_Contact_M(i) = 0;

        case 'D'
            respawnMTB(vec, psi + pi, envelope);
            MT_Growing_D(i) = true;
            MT_GrowthVel_D(i) = Vg;
            MT_Contact_D(i) = 0;
    end
end

function boundaryReached = checkBoundary()
    % Estimate if the pronucleus has impacted the cortex
    global proNucPos R1_max R2_max Prad;

    dist = sqrt((proNucPos(1)/(R1_max - Prad - 0.5))^2 + (proNucPos(2)/(R2_max - Prad - 0.5))^2);
    boundaryReached = dist >= 1;
end

function handleMTGrowth()
    % Handle the growth and shrinkage of microtubules
    global MT_Pos_M MT_Pos_D basePosM basePosD MT_Growing_M MT_Growing_D ...
           MT_GrowthVel_M MT_GrowthVel_D MT_Contact_M MT_Contact_D envelopeM envelopeD ...
           Pr_catastrophe mag_M;

    % Handle MTs for Mother centrosome
    for i = 1:length(MT_Pos_M)
        vecM = MT_Pos_M(i, :) - basePosM;
        mag_M = norm(vecM);

        angleM = atan2(vecM(2), vecM(1));
        if angleM < 0
            angleM = angleM + 2 * pi;
        end

        if mag_M < 0.1
            if MT_Contact_M(i) > 0
                removeContact(angleM);
            end
            respawnMT('M', vecM, i, envelopeM);
        end

        advanceMT(MT_GrowthVel_M(i), vecM, mag_M);
        MT_Pos_M(i, :) = basePosM + vecM;

        if MT_Growing_M(i)
            if rand < Pr_catastrophe && mag_M > 0.4
                MT_Growing_M(i) = false;
                MT_GrowthVel_M(i) = -Vs;
            end
        else
            if rand < Pr_rescue && ~(MT_Contact_M(i) > 0)
                MT_Growing_M(i) = true;
                MT_GrowthVel_M(i) = Vg;
            end
        end

        magScaledM = sqrt((MT_Pos_M(i, 1)/R1_max)^2 + (MT_Pos_M(i, 2)/R2_max)^2);

        if MT_Contact_M(i) > 0
            MT_Contact_M(i) = MT_Contact_M(i) - Tau;
            if MT_Contact_M(i) <= 0
                MT_Contact_M(i) = 0;
                MT_GrowthVel_M(i) = -Vs_c;
                angleM = atan2(MT_Pos_M(i, 2), MT_Pos_M(i, 1));
                if angleM < 0
                    angleM = angleM + 2 * pi;
                end
                removeContact(angleM);
            end
        elseif magScaledM >= 1
            mtContactTest('M', i);
        end
    end

    % Handle MTs for Daughter centrosome
    for i = 1:length(MT_Pos_D)
        vecD = MT_Pos_D(i, :) - basePosD;
        mag_D = norm(vecD);

        angleD = atan2(vecD(2), vecD(1));
        if angleD < 0
            angleD = angleD + 2 * pi;
        end

        if mag_D < 0.1
            if MT_Contact_D(i) > 0
                removeContact(angleD);
            end
            respawnMT('D', vecD, i, envelopeD);
        end

        advanceMT(MT_GrowthVel_D(i), vecD, mag_D);
        MT_Pos_D(i, :) = basePosD + vecD;

        if MT_Growing_D(i)
            if rand < Pr_catastrophe && mag_D > 0.4
                MT_Growing_D(i) = false;
                MT_GrowthVel_D(i) = -Vs;
            end
        else
            if rand < Pr_rescue && ~(MT_Contact_D(i) > 0)
                MT_Growing_D(i) = true;
                MT_GrowthVel_D(i) = Vg;
            end
        end

        magScaledD = sqrt((MT_Pos_D(i, 1)/R1_max)^2 + (MT_Pos_D(i, 2)/R2_max)^2);
%%
        % if MT_Contact_D(i) > 0
        %     MT_Contact_D(i) = MT_Contact_D(i) - Tau;
        %     if MT_Contact_D(i) <= 0
        %         MT_Contact_D(i) = 0;
        %         MT_GrowthVel_D(i) = -Vs_c;
        %         angleD = atan2(MT_Pos_D(i, 2); %, MT_Pos
                
% % Here's the continuation and completion of the MATLAB translation for the `main.cpp` file:
% % 
% % ### 3. `main.m` (Continued)
% % 
% % ```matlab
%%
    % Calculate scaled magnitudes to check contact
    magScaledD = sqrt((MT_Pos_D(i, 1)/R1_max)^2 + (MT_Pos_D(i, 2)/R2_max)^2);

    % Update Contact Indicators
    if MT_Contact_D(i) > 0
        MT_Contact_D(i) = MT_Contact_D(i) - Tau;
        if MT_Contact_D(i) <= 0
            MT_Contact_D(i) = 0;
            MT_GrowthVel_D(i) = -Vs_c;
            angleD = atan2(MT_Pos_D(i, 2), MT_Pos_D(i, 1));
            if angleD < 0
                angleD = angleD + 2 * pi;
            end
            removeContact(angleD);
        end
    elseif magScaledD >= 1
        mtContactTest('D', i);
    end
end

end

function formatted = formatOutput(data, ONLY_COMMA)
    % Format data for output
    if ONLY_COMMA
        separator = ',';
    else
        separator = '|';
    end

    if isvector(data)
        formatted = sprintf(['%f' separator], data);
    elseif ismatrix(data)
        formatted = '';
        [rows, cols] = size(data);
        for r = 1:rows
            rowString = sprintf(['%f' separator], data(r, :));
            formatted = [formatted rowString(1:end-1) '\n']; % Remove trailing separator
        end
    end
end
