%% combine pysch data
% assumes there are pysch_data_1 and psych_data_2

psych_data = [];
for i = 1:numel(psych_data_1) % for each axis
    psych_data{i} = [];
    for j = 1:numel(psych_data_1{i}) % for each condition

        psych_data{i}(j).trial_ids = [psych_data_1{i}(j).trial_ids, psych_data_2{i}(j).trial_ids];
        psych_data{i}(j).bump_dirs = [psych_data_1{i}(j).bump_dirs];
        psych_data{i}(j).bump_correct = psych_data_1{i}(j).bump_correct + psych_data_2{i}(j).bump_correct;
        psych_data{i}(j).bump_total = psych_data_1{i}(j).bump_total + psych_data_2{i}(j).bump_total;

        psych_data{i}(j).psych_curve_data = psych_data{i}(j).bump_correct./psych_data{i}(j).bump_total;
        psych_data{i}(j).psych_curve_data(psych_data{i}(j).bump_dirs > 90) = 1-psych_data{i}(j).psych_curve_data(psych_data{i}(j).bump_dirs > 90);

        % fit psychometric curve
        fit_opts = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [.5 -0.5 .1 100],...
            'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10E-8,'TolX',10E-8,...
            'Weights',psych_data{i}(j).bump_total/sum(psych_data{i}(j).bump_total));
        ft = fittype('a+b*(erf(c*(x-d)))','options',fit_opts);

        psych_fit = [];
        try
            [psych_data{i}(j).psych_fit.fitObj,psych_data{i}(j).psych_fit.gof] = fit(psych_data{i}(j).bump_dirs', psych_data{i}(j).psych_curve_data', ft);
        catch
        end
    end
end
