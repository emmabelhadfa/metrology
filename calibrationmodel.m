function calibrationmodel(TinS_file, Q_file, QinS_file, R_file, RinS_file)
    %% Set up MATLAB %%
    clear;
    format long g;
    close all;
    clc;

    %% Define target points %%
    %because the target is 3 m away and the field of view is 45 degrees by 45
    %degrees, the maximum size of the target is slightly less than 6 m across.
    %To make life simple, we have made the board 4 m across, with the centre
    %LED being defined as (0, 0, 0) in its reference frame

    %note the following nomenclature: t = top, b = bottom, m = middle, r =
    %right, l = left

    board_length = 4;
    tolerance = 0; %depending on the manufacturer. for now, set to 0
    hole_distance = (board_length/2);% + rand(tolerance);

    % ---- IN TARGET COORDINATE FRAME ---- %
    T_tl = [-(hole_distance); (hole_distance);  0; 1];  T_tm = [0; (hole_distance);  0; 1]; T_tr = [(hole_distance); (hole_distance);  0; 1];
    T_ml = [-(hole_distance); 0;                0; 1];  T_mm = [0; 0;                0; 1]; T_mr = [(hole_distance); 0;                0; 1];
    T_bl = [-(hole_distance); -(hole_distance); 0; 1];  T_bm = [0; -(hole_distance); 0; 1]; T_br = [(hole_distance); -(hole_distance); 0; 1];

    T = [T_tl T_tm T_tr T_ml T_mm T_mr T_bl T_bm T_br]; % This is the 3x9 matrix that defines the position of the target points in its own coordinate frame

    % ---- IN SURVEY COORDINATE FRAME ---- %
    T_s = readmatrix(TinS_file); %this is provided as an input

    %% Define test bench points %%
    %In order to fix the test bench in 6 DOF, three measurements are required.
    %Each of these three measurements (test bench points A, B, C) each have
    %their own x, y, z coordinates. The resulting matrices below are therefore
    %3 x 3 matrices


    % ---- IN TEST BENCH COORDINATE FRAME ---- %
    Q = readmatrix(Q_file);

    % ---- IN SURVEY COORDINATE FRAME ---- %
    Q_s = readmatrix(QinS_file); %this is provided as an input

    %% Define Sensor points %%
    %In order to fix the Sensor in 6 DOF, three measurements are required.
    %Each of these three measurements (Sensor points A, B, C) each have
    %their own x, y, z coordinates. The resulting matrices below are therefore
    %3 x 3 matrices

    % ---- IN LWR COORDINATE FRAME ---- %
    R = readmatrix(R_file);

    % ---- IN SURVEY COORDINATE FRAME ---- %
    R_s = readmatrix(RinS_file); %this is provided as an input


    %% Define S to T %%

    % Find the middle of the two point clouds
    Ct = sum(T(1:3,:),2)/size(T,2);
    Cts = sum(T_s(1:3,:),2)/size(T_s,2);

    % Subtract the middles from the point clouds
    Tp = T(1:3,:) - Ct;
    Tsp = T_s(1:3,:) - Cts;

    % Calculate the 'H' matrix and get the single value decomposition
    Ht = Tsp*Tp';
    [Ut,Dt,Vt] = svd(Ht);

    % Ignore the 'D' from the SVD and that is your rotation that best matches
    % the two point clouds together
    Rt = Ut*Vt';

    % This is the translation vector between the two point clouds
    tt = Cts - Rt*Ct;

    % This is the matrix that lines up the 'T' with the 'S' points
    Mt = [Rt tt; 0 0 0 1];

    % Transform the points for the plot
    Smt = Mt*T;

    % Show a figure of the results
    figure;
    hold on;
    set(gcf,"Position",[330 100 940 870]);
    plot3(T(1,:),  T(2,:),  T(3,:),  '*r');
    plot3(T_s(1,:),  T_s(2,:),  T_s(3,:),  'og');
    plot3(Smt(1,:), Smt(2,:), Smt(3,:), '.b');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    legend('Target points', 'Sensor points', 'Target points transformed');

    %% Define S to Q %

    % Find the middle of the two point clouds
    Cq = sum(Q(1:3,:),2)/size(Q,2);
    Cqs = sum(Q_s(1:3,:),2)/size(Q_s,2);

    % Subtract the middles from the point clouds
    Qp = Q(1:3,:) - Cq;
    Qsp = Q_s(1:3,:) - Cqs;

    % Calculate the 'H' matrix and get the single value decomposition
    Hq = Qsp*Qp';
    [Uq,Dq,Vq] = svd(Hq);

    % Ignore the 'D' from the SVD and that is your rotation that best matches
    % the two point clouds together
    Rq = Uq*Vq';

    % This is the translation vector between the two point clouds
    tq = Cqs - Rq*Cq;

    % This is the matrix that lines up the 'Q' with the 'S' points
    Mq = [Rq tq; 0 0 0 1];

    % Transform the points for the plot
    Smq = Mq*Q;

    % Show a figure of the results
    figure;
    hold on;
    set(gcf,"Position",[330 100 940 870]);
    plot3(Q(1,:),  Q(2,:),  Q(3,:),  '*r');
    plot3(Q_s(1,:),  Q_s(2,:),  Q_s(3,:),  'og');
    plot3(Smq(1,:), Smq(2,:), Smq(3,:), '.b');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    legend('Test bench points', 'Total Station points', 'Test bench points transformed');

    %% Define the Transformation from Q to T %%
    Mqt = Mq * inv(Mt);

    %% Define the Points of T in Q %%
    T_inq = Mqt*T; %this gives the x, y, z 

    %% Define S to R %

    % Find the middle of the two point clouds
    Cr = sum(R(1:3,:),2)/size(R,2);
    Crs = sum(R_s(1:3,:),2)/size(R_s,2);

    % Subtract the middles from the point clouds
    Rp = R(1:3,:) - Cr;
    Rsp = R_s(1:3,:) - Crs;

    % Calculate the 'H' matrix and get the single value decomposition
    Hr = Rsp*Rp';
    [Ur,Dr,Vr] = svd(Hr);

    % Ignore the 'D' from the SVD and that is your rotation that best matches
    % the two point clouds together
    Rr = Ur*Vr';

    % This is the translation vector between the two point clouds
    tr = Crs - Rr*Cr;

    % This is the matrix that lines up the 'R' with the 'S' points
    Mr = [Rr tr; 0 0 0 1];

    % Transform the points for the plot
    Smr = Mr*r;

    % Show a figure of the results
    figure;
    hold on;
    set(gcf,"Position",[330 100 940 870]);
    plot3(R(1,:),  R(2,:),  R(3,:),  '*r');
    plot3(R_s(1,:),  R_s(2,:),  R_s(3,:),  'og');
    plot3(Smr(1,:), Smr(2,:), Smr(3,:), '.b');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    legend('Sensor points', 'Total Station points', 'Sensor points transformed');

    %% Define the Transformation from R to T %%
    Mrt = Mr * inv(Mt);

    %% Define the points of T in R %%
    T_inr = Mrt*T; %this gives the x, y, z

    %% Save outputs to files %%
    % The code is meant to save the position of the target in the LWR and test
    % bench frames, respectively. Therefore, it will write this data on two
    % seperate files, which will then be used to determine the calibration
    % uncertainty using uncert_calibrationmodel.m

    % ----- WRITE TARGET IN LWR FILE -----%
    if isfile('target_sensor.csv')
        writematrix(T_inr, 'target_sensor.csv', 'WriteMode', 'append');
    else
        writematrix(T_inr, 'target_sensor.csv');
    end

    % ----- WRITE TARGET IN TEST BENCH FILE -----%
    if isfile('target_testbench.csv')
        writematrix(T_inq, 'target_testbench.csv', 'WriteMode', 'append');
    else
        writematrix(T_inq, 'target_testbench.csv');
    end

end


%% FOR TESTING PURPOSES ONLY %%

% Parameters for a known transform of the Target in the Sensor frame
%rx = 30*(pi/180);
%ry = 6*(pi/180);
%rz = 9*(pi/180);
%tx = 10;
%ty = 5;
%tz = 300;

% Rotation around z
%Rz = [ cos(rz) -sin(rz)  0;
%       sin(rz)  cos(rz)  0;
%       0        0        1];

% Rotation around y
%Ry = [ cos(ry)  0        sin(ry);
%       0        1        0;
%      -sin(ry)  0        cos(ry)];

% Rotation around x
%Rx = [ 1        0        0;
%       0        cos(rx) -sin(rx);
%       0        sin(rx)  cos(rx)];

% Note the order - a different order corresponds to different rx, ry, rz
% values
% R0 = Rz*Ry*Rx;

% Matrix that transforms the target points into the sensor frame
% M0 =  [  R0 [tx; ty; tz];
%               0 0 0 1];

% Transform the points.  When we do this 'for real' we will load these
% points from a file.  This way we can see if M lines up with M0 for
% testing the code.
% S = M0*T;
