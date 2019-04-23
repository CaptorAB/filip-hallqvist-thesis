clear;

max_rows = 2530;
span = 30;

domesticEquityData = readmatrix('../data/omxs301d.csv');
domesticEquityChanges = domesticEquityData(:, 3);
domesticEquityChanges = domesticEquityChanges(1:span:max_rows, :);
domesticEquityModel = fitdist(domesticEquityChanges, 'Normal');

globalEquityData = readmatrix('../data/msciworld1d-filtered.csv');
globalEquityChanges = globalEquityData(:, 3);
globalEquityChanges = globalEquityChanges(1:span:max_rows, :);
globalEquityModel = fitdist(globalEquityChanges, 'Normal');

alternativeChanges = zeros(max_rows, 1);
alternativeChanges = alternativeChanges(1:span:max_rows, :);

creditChanges = zeros(max_rows, 1);
creditChanges = creditChanges(1:span:max_rows, :);

usdSekData = readmatrix('../data/usdsek1d.csv');
usdSekValues = usdSekData(:, 2);
usdSekChanges = diff(usdSekValues)./usdSekValues(1:size(usdSekValues)-1); 
usdSekChanges = usdSekChanges(1:span:max_rows, :);
usdSekModel = fitdist(usdSekChanges, 'Normal');

swapData = readmatrix('../data/swap1d-filtered.csv');

swap2YValues = swapData(:, 3);
swap2YChanges = diff(swap2YValues)./swap2YValues(1:size(swap2YValues)-1); 
swap2YChanges = swap2YChanges(~any(isinf(swap2YChanges), 2), :);
swap2YChanges = swap2YChanges(1:span:max_rows, :);
swap2YModel = fitdist(swap2YChanges, 'Normal');

swap5YValues = swapData(:, 6);
swap5YChanges = diff(swap5YValues)./swap5YValues(1:size(swap5YValues)-1); 
swap5YChanges = swap5YChanges(1:span:max_rows, :);
swap5YModel = fitdist(swap5YChanges, 'Normal');

swap20YValues = swapData(:, 14);
swap20YChanges = diff(swap20YValues)./swap20YValues(1:size(swap20YValues)-1); 
swap20YChanges = swap20YChanges(1:span:max_rows, :);
swap20YModel = fitdist(swap20YChanges, 'Normal');

domesticEquityModel
globalEquityModel
usdSekMode
swap2YModel
swap5YModel
swap20YModel

C = corrcoef([
  domesticEquityChanges';
  globalEquityChanges';
  alternativeChanges';
  swap2YChanges';
  swap5YChanges';
  swap20YChanges';
  creditChanges';
  usdSekChanges'
]');


C