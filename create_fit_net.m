function net = create_fit_net(inputs,targets)
%CREATE_FIT_NET Creates and trains a fitting neural network.
%
%  NET = CREATE_FIT_NET(INPUTS,TARGETS) takes these arguments:
%    INPUTS - RxQ matrix of Q R-element input samples
%    TARGETS - SxQ matrix of Q S-element associated target samples
%  arranged as columns, and returns these results:
%    NET - The trained neural network
%
%  For example, to solve the Simple Fit dataset problem with this function:
%
%    load simplefit_dataset
%    net = create_fit_net(simplefitInputs,simplefitTargets);
%    simplefitOutputs = sim(net,simplefitInputs);
%
%  To reproduce the results you obtained in NFTOOL:
%
%    net = create_fit_net(train_data',test_data');

% Create Network
load train_file.mat
load test_file.mat
inputs=train_data;
targets=test_data;
numHiddenNeurons = 20;  % Adjust as desired
net = newfit(inputs,targets,numHiddenNeurons);
net.divideParam.trainRatio = 70/100;  % Adjust as desired
net.divideParam.valRatio = 15/100;  % Adjust as desired
net.divideParam.testRatio = 15/100;  % Adjust as desired

% Train and Apply Network
[net,tr] = train(net,inputs,targets);
outputs = sim(net,inputs);

% Plot
plotperf(tr)
plotfit(net,inputs,targets)
plotregression(targets,outputs)
