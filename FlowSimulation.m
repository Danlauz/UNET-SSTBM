function [R]=FlowSimulation(nx,ny,Z,BC,cR,Q)

%Source Localisation and Index
x0Q=[ nx/2 ny/2];
LocInj=nx*(x0Q(:,2)-1)+x0Q(:,1);

%% Grid, petrophysics, and fluid objects
% The grid and rock model
G    = computeGeometry(cartGrid([nx ny],[100 100]));
rock = makeRock(G, darcy, 0.3);

b=3; % MOdelling thickness

rock.perm=10.^(Z(1:nx*ny)+log10(3*darcy))*b;

% Fluid properties
fluid = initSimpleADIFluid('phases','W',           ... % Fluid phase: water
                           'mu',  1*centi*poise,   ... % Viscosity
                           'rho', 1000,            ... % Surface density [kg/m^3]
                           'c',   4.4e-10*b,      ... % Fluid compressibility
                           'cR',  cR*b       ... % Rock compressibility
                           );

% Make Reservoir Model
gravity reset off
wModel = WaterModel(G, rock, fluid,'gravity',[0 0 0]);
% Prepare the model for simulation.
wModel = wModel.validateModel();

% Drive mechansims and schedule

% Well: at the midpoint of the south edge
wc = sub2ind(G.cartDims, LocInj);
W = addWell([], G, rock,  wc,     ...
        'Type', 'rate', 'Val', Q, ...
        'Radius', 0.1, 'Name', 'P1','Comp_i',1,'sign',-1);

% Boundary conditions: fixed pressure at E-W sides and no-flow elsewhere
bc  = pside([], G, 'EAST', 9.80638*1000*BC(1), 'sat', 1);
bc  = pside(bc, G, 'WEST', 9.80638*1000*BC(2), 'sat', 1);


% Schedule: describing time intervals and corresponding drive mechanisms
schedule = simpleSchedule(diff(linspace(0,7*day,18)), 'bc', bc, 'W', W);
schedule.step.val=([0.25 0.5 0.75 1 2 5 10 15 20 30 45 60 90 120 180 240 720]-[0 0.25 0.5 0.75 1 2 5 10 15 20 30 45 60 90 120 180 240])*60;

% Reservoir state
state = initResSol(G, 9.80638*1000*(repmat(1:-1/(nx-1):0,1,ny)'*(BC(2)-BC(1))+BC(1)),1);
state.wellSol  = initWellSolAD([], wModel, state);      % No well initially

% Vertical equilibrium
verbose = false;
nonlinear = NonLinearSolver();
stateInit = nonlinear.solveTimestep(state, 10000*day, wModel, 'bc', bc);

% Run simulations
% Simulation pressure
[~, states] = simulateScheduleAD(stateInit, wModel, schedule);

R0=state.pressure*0.000101974;
Rinit=stateInit.pressure*0.000101974;
for k1=1:length(states)
    R(:,k1)=states{k1}.pressure*0.000101974;
end

select=[4 6 8 10 12 14 15 17];
R=[R0 Rinit R(:,select)];


