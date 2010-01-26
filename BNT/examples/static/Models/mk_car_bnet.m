function bnet = mk_car_bnet()
% MK_CAR_BNET Make the car trouble-shooter bayes net.
%
% This network is from p13 of "Troubleshooting under uncertainty", Heckerman, Breese and
% Rommelse, Microsoft Research Tech Report 1994.


BatteryAge = 1;
Battery = 2;
Starter = 3;
Lights = 4;
TurnsOver = 5;
FuelPump = 6;
FuelLine = 7;
FuelSubsys =8;
Fuel = 9;
Spark = 10;
Starts = 11;
Gauge = 12;

n = 12;
dag = zeros(n);
dag(1,2) = 1;
dag(2,[4 5])=1;
dag(3,5) = 1;
dag(6,8) = 1;
dag(7,8) = 1;
dag(8,11) = 1;
dag(9,12) = 1;
dag(10,11) = 1;

arity = 2;
ns = arity*ones(1,n);
bnet = mk_bnet(dag, ns);
for i=1:n
  bnet.CPD{i} = tabular_CPD(bnet, i);
end

  
