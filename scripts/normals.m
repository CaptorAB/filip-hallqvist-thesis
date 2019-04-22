clear;

max_rows = 2530;

domesticEquityData = readmatrix('../data/omxs301d.csv');
domesticEquityChanges = domesticEquityData(:, 3);
domesticEquityModel = fitdist(domesticEquityChanges, 'Normal');

globalEquityData = readmatrix('../data/msciworld1d-filtered.csv');
globalEquityChanges = globalEquityData(:, 3);
globalEquityModel = fitdist(globalEquityChanges, 'Normal');

alternativeChanges = zeros(max_rows, 1);
creditChanges = zeros(max_rows, 1);

usdSekData = readmatrix('../data/usdsek1d.csv');
usdSekValues = usdSekData(:, 2);
usdSekChanges = diff(usdSekValues)./usdSekValues(1:size(usdSekValues)-1); 
usdSekModel = fitdist(usdSekChanges, 'Normal');

swapData = readmatrix('../data/swap1d-filtered.csv');

swap2YValues = swapData(:, 3);
swap2YChanges = diff(swap2YValues)./swap2YValues(1:size(swap2YValues)-1); 
swap2YChanges = swap2YChanges(~any(isinf(swap2YChanges), 2), :);
swap2YModel = fitdist(swap2YChanges, 'Normal');

swap5YValues = swapData(:, 6);
swap5YChanges = diff(swap5YValues)./swap5YValues(1:size(swap5YValues)-1); 
swap5YModel = fitdist(swap5YChanges, 'Normal');

swap20YValues = swapData(:, 14);
swap20YChanges = diff(swap20YValues)./swap20YValues(1:size(swap20YValues)-1); 
swap20YModel = fitdist(swap20YChanges, 'Normal');

swap2YModel
swap5YModel
swap20YModel

C = corrcoef([
  domesticEquityChanges(1:max_rows, 1)';
  globalEquityChanges(1:max_rows, 1)';
  alternativeChanges(1:max_rows, 1)';
  swap2YChanges(1:max_rows, 1)';
  swap5YChanges(1:max_rows, 1)';
  swap20YChanges(1:max_rows, 1)';
  creditChanges(1:max_rows, 1)';
  usdSekChanges(1:max_rows, 1)'
]');


C