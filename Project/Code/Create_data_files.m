
% Clear Workspace
clear

%Load Data
Postmerger_db1b = readtable('postmerger_data_1.csv');
Postmerger_MkDist = readtable('postmerger_data_2.csv');
Premerger_db1b  = readtable('premerger_data_1.csv');
Premerger_MkDist = readtable('premerger_data_2.csv');

%Save files
save("Premerger_db1b.mat", "Premerger_db1b")
save("Premerger_MkDist.mat", "Premerger_MkDist")
save("Postmerger_db1b.mat", "Postmerger_db1b")
save("Postmerger_MkDist.mat", "Postmerger_MkDist")