%%%%%%%%
% Implemented after GSEA --> check the p-value of the top features in asthma dataset
% Make boxplot for the expression log FC values of both diseases
%%%%%%%%

clf; clear;
% Read file
fid = fopen('final_feature_set.txt');
tline = fgetl(fid);
features = {};
count = 1;
while ischar(tline)
    features{count} = tline;
    count = count + 1;
    tline = fgetl(fid);
end
fclose(fid);

% fid = fopen('LC_UP_GENE_SET.tsv');
% C = textscan(fid, '%s %s %s %f %f %f %s', 'HeaderLines', 1);
% fclose(fid);
% features = C{2};

%% asthma
data_asthma = readtable('Asthma_all_genes_p_logfc.csv');

feat_data = data_asthma.ID;
[common_features,idx_common_features] = intersect(feat_data,features,'stable');

new_table = data_asthma(idx_common_features,:);

% log fold change vs label
new_table2 = new_table(new_table.P_Value < 0.05, :);

%% Lung cancer
data_LC = readtable('LC_all_genes_p_logfc.csv');

feat_data2 = data_LC.ID;
[common_features2,idx_common_features2] = intersect(feat_data2,features,'stable');

new_table_LC = data_LC(idx_common_features2,:);
new_table_LC2 = new_table_LC(new_table_LC.P_Value < 0.05, :);

%% Get log FC values in matrix

logFC_A = new_table2.logFC;
logFC_LC = new_table_LC2.logFC;
% subsample LC
n = randperm(size(logFC_A,1));
logFC_LC = logFC_LC(n);

logFC = [logFC_A, logFC_LC];

%% Significance

[h, p] = ttest2(logFC_A,logFC_LC);

%% Box plot
f = figure(1);
hold on;
scatter([1],logFC_A,20,'filled','r');
scatter([2],logFC_LC,20,'filled','b');
boxplot(logFC,'Labels',{'Asthma','Lung Cancer'});
% xlabel({"Asthma","Lung Cancer"});
ylabel("Log Fold Expression Change");
grid on;
text(0.7, 5, "P-value < 0.001", 'Color','k','FontSize',18);
saveas(f,'boxplot_logFC.png');

%% Get top 25 genes --> 22 genes --> couldn't identify gene symbols

new_set_of_features = new_table2.Gene;
new_set_of_features = new_set_of_features(~cellfun('isempty',new_set_of_features));
T = table(new_set_of_features, 'VariableNames', {'Gene_Symbol'} );
writetable(T, 'pruned_features.txt');


