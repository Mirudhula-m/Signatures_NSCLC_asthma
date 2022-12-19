%%%%%%%%%%%%%
% Decision Classification Tree Implementation
%%%%%%%%%%%%%

%% Loading and transforming data

clear
% Load list of genes
genes_classifier = textread('pruned_features.txt','%s');
% Load expression data
asthma = readtable('exprs/asthma/GSE63142_exprs.csv');
LC = readtable('exprs/lc/GSE68793_exprs.csv');
% Load clinical information
asthma_p = readtable('pData/asthma/GSE63142_pData.csv');
lc_p = readtable('pData/lc/GSE68793_pData.csv');
% Load gene annotations
asthma_a = readtable('annot/asthma/GSE63142_annot.csv');
lc_a = readtable('annot/lc/GSE68793_annot.csv');

% capture only pruned feature set

[set2, set2_idx] = intersect(lc_a.GeneSymbol,genes_classifier);
new_genes_classifier = intersect(genes_classifier,set2);
[set1, set1_idx] = intersect(asthma_a.GENE_SYMBOL,new_genes_classifier);

% get ID_ref
new_asthma_a = asthma_a(set1_idx,:);
asthma_IDs = new_asthma_a.ID;

new_lc_a = lc_a(set2_idx,:);
lc_IDs = new_lc_a.ID;

% at this point there are 15 genes in the classifier.
[n, asthma_rows] = intersect(asthma.ID_REF, asthma_IDs);
[n, lc_rows] = intersect(LC.ID_REF, lc_IDs);
asthma2 = asthma(asthma_rows,:);
lc1 = LC(lc_rows,:);

asthma1 = table;
for i = 1:size(asthma2,1)
    asthma1(i,:) = asthma2(find(strcmp(asthma2.ID_REF, asthma_IDs(i))),:);
end
lc2 = table;
for i = 1:size(asthma1,1)
    lc2(i,:) = lc1(find(strcmp(lc1.ID_REF, lc_IDs(i))),:);
end

% asthma1 = asthma1(2:end,:);
% lc2 = lc2(2:end,:);
% Asthma samples for training
[asthma_P,asthma_P_ID] = find(strcmp(table2cell(asthma_p(strcmp(asthma_p.Metadata, 'individual'), :)),'severe asthmatic'));
[lc_P,lc_P_ID] = find(strcmp(table2cell(lc_p(strcmp(lc_p.Metadata, 'smoking history'), :)),'Current reformed smoker for > 15 years'));

% training samples for classification
asthma_train_table = asthma1(:,asthma_P_ID);
lc_train_table = lc2(:,lc_P_ID);

asthma_train = table2array(asthma_train_table); % patient as rows and genes as columns 
lc_train = table2array(lc_train_table);

asthma_z = (asthma_train - mean(asthma_train,2))./(std(asthma_train')');
lc_z = (lc_train - mean(lc_train,2))./(std(lc_train')');


data_exp = [asthma_z,lc_z];
zdata = data_exp;

% Zscore transformation of the data
% zdata = (data_exp - mean(data_exp,2))./(std(data_exp')');

% Discretizing data
threshold = .5;
zdata(zdata<=-1*threshold)=-1;
zdata(zdata>-1*threshold & zdata<=threshold)=0;
zdata(zdata>threshold)=1;

class = ([repmat({'Asthma'},length(asthma_train),1); repmat({'Lung_Cancer'},length(lc_train),1)]);


%% Removing a set of subjects for validation later
idx = randperm(size(zdata, 2));
zdata = zdata(:, idx);
class = class(idx);
val_data = zdata(:,77:end);
class_val = class(77:end);
class_train = class(1:76);
zdata_train = zdata(:,1:76);

%% Running Decision Tree for all possible combinations
combNum = 8;
comb1 = nchoosek(1:size(zdata_train,1),combNum);
for i = 1:size(comb1,1)
    genes = comb1(i,:);
    t = fitctree(zdata_train(genes,:)',class_train,'CrossVal','on','KFold',5);
    dtCVErr(i)  = kfoldLoss(t);
    fprintf("%d of %d\n",i,size(comb1,1))
end

%% Tree Pruning

% sort combinations from lowest to highest error
[IA,IB] = sort(dtCVErr);
% Best combinations
genes = comb1(IB(1),:);
best_set = new_genes_classifier(genes);
best_tree_train = fitctree(zdata_train(genes,:)',class_train);
resubLoss_train = resubLoss(best_tree_train);

best_tree_val = fitctree(val_data(genes,:)',class_val);
resubLoss_val = resubLoss(best_tree_val);
% p_tree = prune(best_tree);
% view(p_tree,'mode','graph');
% view(best_tree,'mode','graph');

clf;
fig = figure;
a = [0.2632,0.1842,0.2237,0.1447,0.2237,0.1316,0.1579,0.1842,0.1447,0.1974,0.1316];
b = [0.2632,0.3158,0.2632,0.2105,0.2105,0.2105,0.2632,0.2632,0.1579,0.2105,0.2632];
r = fit([2:12]',a','poly2');
gg = plot(r,[2:12],a);
set(gg,'Color','b');
hold on;
q = fit([2:12]',b','poly2');
tt = plot(q,[2:12],b);
set(tt,'Color','r');
legend("Train loss data","Train loss fit curve","Test loss data","Test loss fit curve");
grid on;
ylabel("Loss");
xlabel("Number of features");
saveas(fig,"Loss_curve_Dtrees.png");

