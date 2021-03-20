clc
tic
load monkeydata0
tests = testFunction_for_Classifier('Team Cassius');
% tests = testFunction_for_Classifier_2('Team Cassius');
% tests = testFunction_for_Classifier_3('Team Cassius');
% tests = testFunction_for_Classifier_4('Team Cassius');
% tests = testFunction_for_Classifier_5('Team Cassius');
elapsed_time = toc;
fprintf('Total elapsed time is %f', elapsed_time);
fprintf('Time elapsed per test is %f', elapsed_time/tests);

