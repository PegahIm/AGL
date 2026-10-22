(* ::Package:: *)

(*  Bigram frequency setup*)

(*Training items,anchored at both beginning ("b") and end ("e")*)
strings1 = {"bXXVTe", "bXXVXJJe", "bVXJJe", "bXVJTVJe", "bXXVXJe",
   "bXVXe", "bXXXVTe", "bVJe", "bXVXJJJe", "bVJTVTVe", "bVJTVXe",
   "bXXXVTVe", "bXXVJe", "bXVXJJe", "bVTe", "bVJTVXJe", "bXXXVXe",
   "bVJTXVJe", "bXVXJe", "bXXXXVXe", "bXVTe", "bVJTVJe", "bVXJJJJe"};

strings2 = StringDrop[#, -1] & /@ strings1;
bigramsFunc[strings_] :=
  Map[StringJoin @@@ Partition[Characters[#], 2, 1] &, strings];
bigramCountsFunc[bigrams_] := Tally[Flatten[bigrams]];
bigrams = bigramsFunc[strings2];
bigramCounts = bigramCountsFunc[bigrams];
(*Group the bigrams by their first letter*)
groupedBigrams = GroupBy[bigramCounts, StringTake[First[#], 1] &];
(*Calculate conditional bigram frequencies*)
bigramFrequencies =
  KeyValueMap[
   Function[{firstLetter, bigrams},
    Map[Function[{bigram}, {First[bigram],
       Last[bigram]/Total[Last /@ bigrams]}], bigrams]],
   groupedBigrams];
(*Add terminal bigrams*)
additionalBigrams = {{"Xe", 4/23}, {"Te", 4/23}, {"Ve", 2/23}, {"Je",
    13/23}};
bigramFrequencies = Append[bigramFrequencies, additionalBigrams];

allBigrams = {"bX", "bV", "bJ", "bT", "XX", "XV", "XJ", "VX", "VT",
   "VJ", "TV", "TX", "JJ", "JT", "JV", "Xe", "Ve", "Je", "Te", "JX",
   "XT", "TJ", "VV"};

flatBigramsList = Flatten[bigramFrequencies, 1];
missingBigrams = Complement[allBigrams, flatBigramsList[[All, 1]]];
flatBigramsList = Join[flatBigramsList, Thread[{missingBigrams, 0}]];


(* Model formulas *)

calculateCoefficientsQuantum[freq_?NumericQ, th1_?NumericQ] :=
  Module[{alpha, beta, alphaPrime, betaPrime},
   alpha = Exp[I th1] Sqrt[freq];
   beta = -Sqrt[1 - freq];        
   alphaPrime = Sqrt[1 - freq];   
   betaPrime = Exp[-I th1] Sqrt[freq];
   {{alpha, alphaPrime}, {beta, betaPrime}}];

buildBigramCalculationsQuantum[th1_?NumericQ] :=
  Module[{assoc},
   assoc = Association[
     Map[Function[{bg},
       If[Length[bg] == 2,
        bg[[1]] -> calculateCoefficientsQuantum[N[bg[[2]]], th1],
        bg -> {{0, 1}, {-1, 0}}]], flatBigramsList]];
   Do[If[! KeyExistsQ[assoc, bg], assoc[bg] = {{0, 1}, {-1, 0}}], {bg,
      allBigrams}];
   assoc];

measurementOperatorQuantum[x_] := {{Sqrt[1 - x], 0}, {0, Sqrt[x]}};

calculateQuantumProbability[string_, bigramCalcs_, x_?NumericQ] :=
  Module[{bigramsInStr, firstBigram, stateVector, conjugateStateVector,
    M, U, result},
   bigramsInStr = StringJoin /@ Partition[Characters[string], 2, 1];
   firstBigram = First[bigramsInStr];
   stateVector = {{bigramCalcs[firstBigram][[1, 1]]},
     {bigramCalcs[firstBigram][[2, 1]]}};
   M = measurementOperatorQuantum[x];
   conjugateStateVector = ConjugateTranspose[stateVector];
   Do[U = bigramCalcs[bg];
    conjugateStateVector = conjugateStateVector . ConjugateTranspose[U],
    {bg, Rest[bigramsInStr]}];
   result = conjugateStateVector . M;
   stateVector = {{bigramCalcs[firstBigram][[1, 1]]},
     {bigramCalcs[firstBigram][[2, 1]]}}; (* reset -- consumed above *)
   Do[U = bigramCalcs[bg];
    stateVector = U . stateVector, {bg, Rest[bigramsInStr]}];
   Re[(result . M . stateVector)[[1, 1]]]];

bigramCalculationsClassicalBase =
  Association[
   Map[Function[{bg},
     If[Length[bg] == 2,
      bg[[1]] -> {{N[bg[[2]]], 1 - N[bg[[2]]]}, {1 - N[bg[[2]]],
         N[bg[[2]]]}}, bg -> {{0, 1}, {1, 0}}]], flatBigramsList]];
Do[If[! KeyExistsQ[bigramCalculationsClassicalBase, bg],
   bigramCalculationsClassicalBase[bg] = {{0, 1}, {1, 0}}], {bg, allBigrams}];

bigramCalculationsDelta[delta_?NumericQ] :=
  Association[
   Table[With[{p = bigramCalculationsClassicalBase[bg][[1, 1]]},
     With[{c2top = Clip[1 - p - delta, {0, 1}]},
      bg -> {{p, c2top}, {1 - p, 1 - c2top}}]], {bg, allBigrams}]];

measurementOperatorClassical[x_] := {{1 - x, 0}, {0, x}};

calculateClassicalProbabilityDelta[string_, bigramCalcs_, x_?NumericQ] :=
  Module[{bigramsInStr, firstBigram, stateVector, M, L, U, result},
   bigramsInStr = StringJoin /@ Partition[Characters[string], 2, 1];
   firstBigram = First[bigramsInStr];
   stateVector = {bigramCalcs[firstBigram][[1, 1]],
     bigramCalcs[firstBigram][[2, 1]]};
   Do[U = bigramCalcs[bg];
    stateVector = U . stateVector, {bg, Rest[bigramsInStr]}];
   M = measurementOperatorClassical[x];
   L = {1, 1};
   result = L . M;
   Chop[result . stateVector]];


(* Forward prediction and synthetic data generation *)

quantumPredictions[items_, th1_?NumericQ, x_?NumericQ] :=
  Module[{bc = buildBigramCalculationsQuantum[th1]},
   calculateQuantumProbability[#, bc, x] & /@ items];

classicalPredictions[items_, delta_?NumericQ, x_?NumericQ] :=
  Module[{bc = bigramCalculationsDelta[delta]},
   Flatten[calculateClassicalProbabilityDelta[#, bc, x] & /@ items]];

simulateData[trueProbs_, n_] :=
  RandomVariate[BinomialDistribution[n, #]] & /@
   Clip[trueProbs, {0.0001, 0.9999}];


(* Refitting functions *)

fitBayesianSim[yData_, testItems_, n_, maxIters_: 200] :=
  Module[{nItems = Length[testItems], logLC, opt, xOptC, deltaOptC,
    logLCVal, predictedC, G2C, BICC},
   logLC[xv_?NumericQ, dv_?NumericQ] :=
    Module[{bc, modelVals},
     bc = bigramCalculationsDelta[dv];
     modelVals = calculateClassicalProbabilityDelta[#, bc, xv] & /@ testItems;
     If[! AllTrue[modelVals, 0 < # < 1 &], Return[-Infinity]];
     Sum[Log[Binomial[n, yData[[i]]]] +
       yData[[i]]*Log[modelVals[[i]]] +
       (n - yData[[i]])*Log[1 - modelVals[[i]]], {i, nItems}]];
   opt = Quiet@NMaximize[{logLC[xv, dv], 0 <= xv < 0.5 && -0.5 <= dv <= 0.5},
      {xv, dv}, MaxIterations -> maxIters,
      Method -> {"DifferentialEvolution", "PostProcess" -> False}];
   logLCVal = opt[[1]];
   {xOptC, deltaOptC} = {xv, dv} /. opt[[2]];
   predictedC = classicalPredictions[testItems, deltaOptC, xOptC];
   G2C = 2*Sum[
      yData[[i]]*Log[yData[[i]]/(n*predictedC[[i]])] +
       (n - yData[[i]])*Log[(n - yData[[i]])/(n*(1 - predictedC[[i]]))],
      {i, nItems}];
   BICC = G2C + 2*Log[nItems]; (* k=2: x, delta *)
   <|"x" -> xOptC, "delta" -> deltaOptC, "logL" -> logLCVal, "G2" -> G2C,
     "BIC" -> BICC|>];

fitQuantumSim[yData_, testItems_, n_, maxIters_: 200] :=
  Module[{nItems = Length[testItems], logLQ, opt, theta1OptQ, xOptQ,
    logLQVal, predictedQ, G2Q, BICQ},
   logLQ[th_?NumericQ, xv_?NumericQ] :=
    Module[{modelVals = quantumPredictions[testItems, th, xv]},
     Sum[Log[Binomial[n, yData[[i]]]] +
       yData[[i]]*Re[Log[modelVals[[i]]]] +
       (n - yData[[i]])*Re[Log[1 - modelVals[[i]]]], {i, nItems}]];
   opt = Quiet@NMaximize[{logLQ[th, xv], 0 <= th < 2 Pi && 0 <= xv < 0.5},
      {th, xv}, MaxIterations -> maxIters,
      Method -> {"DifferentialEvolution", "PostProcess" -> False}];
   logLQVal = opt[[1]];
   {theta1OptQ, xOptQ} = {th, xv} /. opt[[2]];
   predictedQ = quantumPredictions[testItems, theta1OptQ, xOptQ];
   G2Q = Re[2*Sum[
       yData[[i]]*Log[yData[[i]]/(n*predictedQ[[i]])] +
        (n - yData[[i]])*Log[(n - yData[[i]])/(n*(1 - predictedQ[[i]]))],
       {i, nItems}]];
   BICQ = G2Q + 2*Log[nItems]; 
   <|"theta1" -> theta1OptQ, "x" -> xOptQ, "logL" -> logLQVal, "G2" -> G2Q,
     "BIC" -> BICQ|>];


(* Test items *)

testItems = {
   "bXXVXJe", "bXVTVJJe", "bVXJe", "bXXVTVe", "bXVJTVXe", "bXXVTVJe",
   "bXXVXe", "bXVTVJe", "bVJTXVTe", "bVJTXVXe", "bVTVJe", "bVJTVTe",
   "bVTVe", "bXVTVe", "bVTVJJe", "bVXe", "bXXJJe", "bXXVe", "bXVXVe",
   "bXVXVJe", "bXXVJJJe", "bXJJe", "bVXVJe", "bXXVVJJe", "bXVXTe",
   "bVXJTJe", "bVXJJVe", "bXXTXe", "bTVJe", "bVXJJXe", "bVJJXVTe",
   "bVJTVe"};


(* Generating parameter grid*)

n = 30; 
nReps = 200; 
maxIters = 200; 

quantumGeneratingPoints = {
   {5.129, 0.315},
   {5.178, 0.306},    
   {5.166, 0.278}    
   };
bayesianGeneratingPoints = {
   {-0.163, 0.0},
   {-0.174, 0.0},
   {-0.171, 0.0}
   };


(* Recovery loop *)

runRecovery[testItemsArg_, conditionLabel_] :=
  Module[{results = {}, trueProbsQ, trueProbsC, yData, fitB, fitQ,
    winner, tStart},

   (* Generate from QUANTUM, refit both *)
   Do[
    trueProbsQ = quantumPredictions[testItemsArg, pt[[1]], pt[[2]]];
    Print["[", conditionLabel, "] Quantum-generated, theta1=", pt[[1]],
      ", x=", pt[[2]], " -- ", nReps, " reps"];
    tStart = AbsoluteTime[];
    Do[
     yData = simulateData[trueProbsQ, n];
     fitB = fitBayesianSim[yData, testItemsArg, n, maxIters];
     fitQ = fitQuantumSim[yData, testItemsArg, n, maxIters];
     winner = If[fitB["BIC"] < fitQ["BIC"], "Bayesian", "Quantum"];
     AppendTo[results,
      <|"Condition" -> conditionLabel, "GeneratingModel" -> "Quantum",
        "TrueTheta1" -> pt[[1]], "TrueX" -> pt[[2]], "Rep" -> rep,
        "FitBayesianX" -> fitB["x"], "FitBayesianDelta" -> fitB["delta"],
        "FitQuantumTheta1" -> fitQ["theta1"], "FitQuantumX" -> fitQ["x"],
        "BayesianBIC" -> fitB["BIC"], "QuantumBIC" -> fitQ["BIC"],
        "Winner" -> winner|>];
     If[Mod[rep, 25] == 0,
      Module[{elapsed = AbsoluteTime[] - tStart, perRep},
       perRep = elapsed/rep;
       Print["[", conditionLabel, "]   rep ", rep, "/", nReps,
        " -- ", NumberForm[perRep, 3], " sec/rep, est. remaining: ",
        NumberForm[perRep*(nReps - rep)/60., 3], " min"]]],
     {rep, 1, nReps}],
    {pt, quantumGeneratingPoints}];

   (* Generate from BAYESIAN, refit both *)
   Do[
    trueProbsC = classicalPredictions[testItemsArg, pt[[1]], pt[[2]]];
    Print["[", conditionLabel, "] Bayesian-generated, delta=", pt[[1]],
      ", x=", pt[[2]], " -- ", nReps, " reps"];
    tStart = AbsoluteTime[];
    Do[
     yData = simulateData[trueProbsC, n];
     fitB = fitBayesianSim[yData, testItemsArg, n, maxIters];
     fitQ = fitQuantumSim[yData, testItemsArg, n, maxIters];
     winner = If[fitB["BIC"] < fitQ["BIC"], "Bayesian", "Quantum"];
     AppendTo[results,
      <|"Condition" -> conditionLabel, "GeneratingModel" -> "Bayesian",
        "TrueDelta" -> pt[[1]], "TrueX" -> pt[[2]], "Rep" -> rep,
        "FitBayesianX" -> fitB["x"], "FitBayesianDelta" -> fitB["delta"],
        "FitQuantumTheta1" -> fitQ["theta1"], "FitQuantumX" -> fitQ["x"],
        "BayesianBIC" -> fitB["BIC"], "QuantumBIC" -> fitQ["BIC"],
        "Winner" -> winner|>];
     If[Mod[rep, 25] == 0,
      Module[{elapsed = AbsoluteTime[] - tStart, perRep},
       perRep = elapsed/rep;
       Print["[", conditionLabel, "]   rep ", rep, "/", nReps,
        " -- ", NumberForm[perRep, 3], " sec/rep, est. remaining: ",
        NumberForm[perRep*(nReps - rep)/60., 3], " min"]]],
     {rep, 1, nReps}],
    {pt, bayesianGeneratingPoints}];

   results];
   
   resultsRecovery = runRecovery[testItems, "AGL-Exp2"]; 


(*Summary tables *)

summarizeIdentification[results_] :=
  Module[{byPoint},
   byPoint =
    GatherBy[results,
     If[#["GeneratingModel"] == "Quantum", {"Quantum", #["TrueTheta1"], #["TrueX"]},
       {"Bayesian", #["TrueDelta"], #["TrueX"]}] &];
   Table[
    Module[{genModel = sub[[1]]["GeneratingModel"], correct, total},
     correct = Count[sub, r_ /; r["Winner"] == genModel];
     total = Length[sub];
     <|"Condition" -> sub[[1]]["Condition"], "GeneratingModel" -> genModel,
       "TrueParams" ->
        If[genModel == "Quantum", {sub[[1]]["TrueTheta1"], sub[[1]]["TrueX"]},
          {sub[[1]]["TrueDelta"], sub[[1]]["TrueX"]}],
       "N" -> total, "PercentCorrectlyRecovered" -> N[100*correct/total]|>],
    {sub, byPoint}]];

summarizeParameterRecovery[results_] :=
  Module[{quantumRows, bayesianRows, byQ, byB},
   quantumRows = Select[results, #["GeneratingModel"] == "Quantum" &];
   bayesianRows = Select[results, #["GeneratingModel"] == "Bayesian" &];
   byQ = GatherBy[quantumRows, {#["TrueTheta1"], #["TrueX"]} &];
   byB = GatherBy[bayesianRows, {#["TrueDelta"], #["TrueX"]} &];
   {
    Table[
     <|"Condition" -> sub[[1]]["Condition"], "GeneratingModel" -> "Quantum",
       "TrueTheta1" -> sub[[1]]["TrueTheta1"], "TrueX" -> sub[[1]]["TrueX"],
       "MeanFitTheta1" -> Mean[#["FitQuantumTheta1"] & /@ sub],
       "RMSETheta1" ->
        Sqrt[Mean[(#["FitQuantumTheta1"] - sub[[1]]["TrueTheta1"])^2 & /@ sub]],
       "MeanFitX" -> Mean[#["FitQuantumX"] & /@ sub],
       "RMSEX" -> Sqrt[Mean[(#["FitQuantumX"] - sub[[1]]["TrueX"])^2 & /@ sub]]|>,
     {sub, byQ}],
    Table[
     <|"Condition" -> sub[[1]]["Condition"], "GeneratingModel" -> "Bayesian",
       "TrueDelta" -> sub[[1]]["TrueDelta"], "TrueX" -> sub[[1]]["TrueX"],
       "MeanFitDelta" -> Mean[#["FitBayesianDelta"] & /@ sub],
       "RMSEDelta" ->
        Sqrt[Mean[(#["FitBayesianDelta"] - sub[[1]]["TrueDelta"])^2 & /@ sub]],
       "MeanFitX" -> Mean[#["FitBayesianX"] & /@ sub],
       "RMSEX" -> Sqrt[Mean[(#["FitBayesianX"] - sub[[1]]["TrueX"])^2 & /@ sub]]|>,
     {sub, byB}]
    }];



(*idTable = summarizeIdentification[resultsRecovery];
{paramTableQ, paramTableB} = summarizeParameterRecovery[resultsRecovery];
Export[FileNameJoin[{NotebookDirectory[], "AGL-Exp2-ModelRecovery-IdentificationTable.csv"}], idTable];
Export[FileNameJoin[{NotebookDirectory[], "AGL-Exp2-ModelRecovery-QuantumParamRecovery.csv"}], paramTableQ];
Export[FileNameJoin[{NotebookDirectory[], "AGL-Exp2-ModelRecovery-BayesianParamRecovery.csv"}], paramTableB];
Print["Exported to: ", NotebookDirectory[]];*)
