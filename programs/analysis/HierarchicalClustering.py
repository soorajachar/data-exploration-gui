#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import shutil
from hal import HAL
import pickle
import os
import re
import copy
from sklearn.preprocessing import MinMaxScaler


class UHAL():

    def __init__(self, data, n_tsne=40000, trn_size=0.8, clf_type='svm', name='default', phenotype_complexity = 2,arcsinh=True, scaling=True):

        self.data = data  # training and testing data
        self.model = None  # HAL model
        self.n_tsne = n_tsne  # maximum number of data points to be used for tSNE transformation and model training
        self.trn_size = trn_size  # fraction of data used for training
        self.clf_type = clf_type  # classifier type (parameter for HAL-x)
        self.phenotypes = None  # stores binary phenotypes
        self.name = name  # name of run
        self.arcsinh = arcsinh  # whether to apply arcsinh transformation
        self.scaling = scaling  # whether to apply standard (z-score) scaling (STRONGLY recommended)
        self.scaler_obj = None  # scaling mechanism (used everywhere to prevent discrepancies)

        self.phenotype_complexity = [] #### Breakdown of digitized phenotypes
        for i in range(1,phenotype_complexity+1):
            self.phenotype_complexity.append(i/phenotype_complexity)
        
        

    def defineModel(self, markers, score, output_dir=None, retrain=False):
        """ This method creates a training set and runs the clustering algorithm

        markers: a list of strings corresponding to columns in the input data frames
        score: score to use in new model
        arcsinh: whether to apply the arcsinh transformation prior to clustering (this is usually a good idea)
        scaling: z-score scaling (this is considered best practice)
        output_dir (default is None): output directory, None indicates the current working directory
        retrain (default is False): specifies whether to retrain with existing model or create a new one
        """
        assert (not retrain) or (self.model is not None)  # exception if retraining without existing model

        if output_dir is None:
            output_dir = os.getcwd()  # set to current directory if not specified

        data = copy.copy(self.data)  # dictionary of samples (file names are keys)

        n_cells = np.floor(self.n_tsne / len(set(data.index.get_level_values(0))))  # number of data points selected from each sample
        samples = list(set(data.index.get_level_values(0)))  # sample/file names
        origins = []  # list tracking which sample individual points came from

        for ii, sample in enumerate(list(set(data.index.get_level_values(0)))):  # iterate through samples
            sample_data = data[data.index.get_level_values(0) == sample]
            sample_size = int(np.min([n_cells, sample_data.shape[0]]))  # number of points to select per sample
            random_choice = np.random.choice(sample_data.shape[0], sample_size)
            origins.extend([sample]*len(random_choice))  # note where data points came from
            if ii == 0:
                data_samples = sample_data[markers].iloc[random_choice, :].values  # start list of data points
            else:
                data_samples = np.concatenate([data_samples, sample_data[markers].iloc[random_choice, :].values])
        '''
        for i, current_marker in enumerate(markers):
            print(current_marker)
            print(stats.entropy(np.arcsinh(data_samples[:, i])))
            plt.hist(np.arcsinh(data_samples[:, i]))
            plt.show()
        '''
        # determine whether the current experiment has been processed (with any CV score)
        redundant = False
        for file_name in os.listdir(output_dir+'/serialized'):
            match = re.fullmatch('model_0.*%s\.pkl' % self.name, file_name)
            if match is not None:
                redundant = True
                break
        # create new model if not retraining
        model_file = 'model_' + self.name + '.pkl'
        scaler_file = 'scaler_' + self.name + '.pkl'
        label_file = 'Labels_tSNE_' + str(score) + self.name + '.pkl'
        if (label_file in os.listdir(output_dir+'/serialized')) and not retrain:
            # re-run experiment with same CV score
            model = pickle.load(open(output_dir + '/serialized/' + model_file, 'rb'))
            self.scaler_obj = pickle.load(open(output_dir + '/serialized/' + scaler_file, 'rb'))
            tsne_frame = pickle.load(open(output_dir + '/serialized/' + label_file, 'rb'))
            labels_tSNE = tsne_frame['clusters']
            data_samples = tsne_frame.loc[:, markers].values
            output = tsne_frame
        else:
            if redundant and not retrain:
                # re-run experiment with different CV score
                model = pickle.load(open(output_dir + '/serialized/' + model_file, 'rb'))
                data_samples = pickle.load(open(output_dir + '/serialized/tSNE_subset_' + self.name + '.pkl', 'rb'))
                self.scaler_obj = pickle.load(open(output_dir + '/serialized/' + scaler_file, 'rb'))
            else:
                # create HAL object and fit model to data (using only training data)
                try:
                    shutil.rmtree('./info_hal')  # remove old info_hal folder
                except FileNotFoundError:
                    pass
                model = HAL(clf_type=self.clf_type, outlier_ratio=0.1, late_exag=900, alpha_late=2.0, n_cluster_init=150,
                            warm_start=True)
            # apply arcsinh transformation (enabled by default)
            if self.arcsinh:
                transformed_samples = np.arcsinh(data_samples)
            else:
                transformed_samples = data_samples
            # apply standard (z-score) scaling (enabled by default)
            if self.scaling:
                if self.scaler_obj is None:
                    self.scaler_obj = MinMaxScaler()
                    scaled_data = self.scaler_obj.fit_transform(transformed_samples)
                else:
                    scaled_data = self.scaler_obj.transform(transformed_samples)
            else:
                scaled_data = transformed_samples  # do not use this option without a good reason!
            model.fit(scaled_data)
            pickle.dump(model, open(output_dir+'/serialized/'+model_file, 'wb'))
            pickle.dump(self.scaler_obj, open(output_dir+'/serialized/'+scaler_file, 'wb'))
            # create a frame with the clusters and samples for each data point
            labels_tSNE = model.predict(scaled_data, cv=score)
            output = pd.DataFrame(data_samples)
            output.columns = markers
            output["clusters"] = labels_tSNE
            output["origin"] = origins
            output = self.addTsne(output)
            output.to_csv(output_dir + '/Labels_tSNE_' + str(score) + self.name + '.csv')
            pickle.dump(data_samples, open(output_dir + '/serialized/tSNE_subset_' + self.name + '.pkl', "wb"))
            pickle.dump(output, open(output_dir + '/serialized/Labels_tSNE_' + str(score) + self.name + '.pkl', "wb"))


        self.model = model
        labels_only = np.array(labels_tSNE)

        return labels_only, output  # do not return samples of origin with labels

    def run(self, all_markers, score, unused=[], output_dir=None, retrain=False):
        """ runs the clustering procedure on the data and puts the results in the specified folder

        markers: a list of strings corresponding to columns in the input data frames
        score: score to use in new model
        output_dir (default is None): location of output data (None indicates the current directory)
        retrain (default is False): specifies whether to retrain with existing model or create a new one
        """
        if output_dir is None:
            output_dir = os.getcwd()
        else:
            try:
                os.mkdir(output_dir)
            except FileExistsError:
                pass
        # create folder for serialized files
        try:
            os.mkdir(output_dir+'/serialized')
        except FileExistsError:
            pass
        # train clustering model
        markers = [item for item in all_markers if item not in unused]
        labels_tSNE, tSNE_data = self.defineModel(markers, score, output_dir=output_dir, retrain=retrain)
        # apply model to whole data set, then write results to files
        cluster_labels, all_labels, phenotypes, binary_phens = self.fitModel(markers, score, output_dir=output_dir,
                                                                             unused=unused)
        cf, ce = self.outputData(cluster_labels, all_markers, score=score, output_dir=output_dir)

        return all_labels, {'labels': all_labels, 'frequencies': cf, 'expression': ce, 'phenotypes': phenotypes,
                'binary_phenotypes': binary_phens, 'tsne_subset': tSNE_data}

    def fitModel(self, markers, score, unused=[], output_dir=None):
        """ uses a pre-trained HAL model to assign all data points to clusters
        data_samples: Numpy array of training data (obtained from defineModel)
        labels_tSNE: cluster labels corresponding to data_samples
        """

        all_markers = self.data.columns
        approved_markers = [item for item in all_markers if item in markers and (item != "Labels")]  # exclude Labels

        # label all data points
        '''
        for s in samples:  # iterate through samples and label all points
            relevant_data = self.data[file][approved_markers].values  # select current sample (and only certain markers)
            if self.arcsinh:
                relevant_data = np.arcsinh(relevant_data)
            if self.scaling:
                relevant_data = self.scaler_obj.transform(relevant_data)
            self.data[file]['Labels'] = self.model.predict(relevant_data, cv=score)  # apply clustering model to whole frame
            labels.extend(self.data[file]['Labels'])
        '''
        test_data = self.data[approved_markers].values
        if self.arcsinh:
            test_data = np.arcsinh(test_data)
        if self.scaling:
            test_data = self.scaler_obj.transform(test_data)
        labels = self.model.predict(test_data, cv=score)
        self.data.loc[:, 'Labels'] = labels
        unique_cluster_labels = list(set(labels))
        Phenotype_define = []

        # get phenotypes
        approved_markers.extend(unused)
        for ii, cluster in enumerate(unique_cluster_labels):  # iterate through clusters
            matches = self.data['Labels'] == cluster
            scaled_data = np.arcsinh(self.data.loc[matches, approved_markers])
            marker_means = np.sinh(np.mean(scaled_data, axis=0))
            Phenotype_define.append(marker_means)  # find marker expression in each cluster
        phenotype_frame = pd.DataFrame(np.array(Phenotype_define), columns=approved_markers, index=unique_cluster_labels)

        # create binary phenotypes (positive/negative)
        binary_phenotypes,list_phenotypes = self.numericPhenotype(phenotype_frame,approved_markers)
        self.phenotypes = list_phenotypes
        # pickle and write data
        pickle.dump(binary_phenotypes,
                    open(output_dir + '/serialized/' + 'BinaryPhenotype_' + self.name + "_" + str(score) + '.pkl', "wb"))
        pickle.dump(phenotype_frame,
                    open(output_dir + '/serialized/' + 'Phenotype_' + self.name + "_" + str(score) + '.pkl', "wb"))
        pickle.dump(self.data, open(
            output_dir + '/serialized/' + 'AllData_training_withLabels_' + self.name + "_" + str(score) + '.pkl', "wb"))

        # write CSV
        phenotype_frame.to_csv(output_dir + '/Phenotype_' + self.name + "_" + str(score) + '.csv')
        binary_phenotypes.to_csv(output_dir + '/BinaryPhenotype_' + self.name + "_" + str(score) + '.csv')
        return unique_cluster_labels, labels, phenotype_frame, binary_phenotypes

    def outputData(self, cluster_labels, markers, score, output_dir=None):
        """Write data to output files

        cluster_labels: list of cluster labels (integers)
        markers: all markers being used in the current analysis
        output_dir (default is None): output directory, None indicates the current working directory
        """

        data = copy.copy(self.data)
        samples = list(set(self.data.index.get_level_values(0)))
        popMean = np.zeros((len(samples), len(cluster_labels)*len(markers)))
        populationFrequencies = np.zeros((len(samples), len(cluster_labels)))

        # create list of vetoed markers
        all_markers = data.columns
        vetoed_markers = [item for item in all_markers if item not in markers and (item != "Labels")]  # exclude Labels
        popMean_veto = np.zeros((len(samples), len(cluster_labels)*len(vetoed_markers)))

        # create column names
        cols = []
        for cluster in cluster_labels:
            for mm in markers:
                cols.append(mm + " dual count for Cluster #{} (Mean)".format(cluster))
        cols_veto = []
        for cluster in cluster_labels:
            for vm in vetoed_markers:
                cols_veto.append(vm + " dual count for Cluster #{} (Mean)".format(cluster))

        for file_index, file_num in enumerate(samples):  # iterate through samples
            sample_data = self.data.iloc[self.data.index.get_level_values(0)==file_num, :]
            for ii, cluster in enumerate(cluster_labels):  # iterate through clusters
                labels = sample_data['Labels']
                if cluster not in list(labels):
                    continue  # skip blank phenotypes
                # determine percentage of sample belonging to a given cluster
                populationFrequencies[file_index, ii] = float(np.sum(sample_data['Labels'].values == cluster) / sample_data.shape[0] * 100)
                for jj, mm in enumerate(markers):
                    # record average expression per marker per cluster per sample
                    if self.arcsinh:
                        relevant_data = np.arcsinh(sample_data[sample_data['Labels'] == cluster][mm].values)
                    else:
                        relevant_data = sample_data[sample_data['Labels'] == cluster][mm].values
                    popMean[file_index, (ii*len(markers)) + jj] = np.mean(relevant_data)
                for jj, vm in enumerate(vetoed_markers):
                    if self.arcsinh:
                        relevant_data = np.arcsinh(sample_data[sample_data['Labels'] == cluster][vm].values)
                    else:
                        relevant_data = sample_data[sample_data['Labels'] == cluster][vm].values
                    popMean_veto[file_index, (ii*len(vetoed_markers)) + jj] = np.mean(relevant_data)

        ##### THIS DATA CAN BE USED TO EVALUATE AN ENTIRE MODEL BASED ON AUC
        clusterFrequencies = pd.DataFrame(populationFrequencies, columns=cluster_labels, index=samples)
        clusterExpression = pd.DataFrame(popMean, columns=cols, index=samples)
        clusterExpression2 = pd.DataFrame(popMean_veto, columns=cols_veto, index=samples)

        clusterFrequencies = self.add_phenotypes(clusterFrequencies.copy())
        clusterExpression = self.add_phenotypes(clusterExpression.copy())
        
        
        # serialize and write to files
        if output_dir is None:
            output_dir = os.getcwd()
        pickle.dump(clusterFrequencies, open('%s/serialized/ClusterFrequencies_%s_%s.pkl' %
                                             (output_dir, self.name, str(score)), "wb"))
        pickle.dump(clusterExpression, open('%s/serialized/ClusterDualCountMean_%s_%s.pkl' %
                                            (output_dir, self.name, str(score)), "wb"))
        pickle.dump(clusterExpression2, open('%s/serialized/vetoedMarkerClusterExpression_%s_%s.pkl' %
                                             (output_dir, self.name, str(score)), "wb"))

        # write to CSV files
        clusterFrequencies.to_csv('%s/ClusterFrequencies_%s_%s.csv' % (output_dir, self.name, str(score)))
        clusterExpression.to_csv('%s/ClusterDualCountMean_%s_%s.csv' % (output_dir, self.name, str(score)))
        clusterExpression2.to_csv('%s/vetoedMarkerClusterExpression_%s_%s.csv' % (output_dir, self.name, str(score)))

        return clusterFrequencies, clusterExpression


    def numericPhenotype(self,df,markers,transpose = False):
        """Define whether marker expression is positive or negative relative to median level"""
        phenotypes = []
        
        if transpose:
            df = df.copy().T
        
        for ind in list(df.index):
            phen = []
            for col in markers:
                phen.append(self.assignRegion(df.loc[ind][col],df[col]))
            phenotypes.append(phen)
                
        return pd.DataFrame(np.array(phenotypes), columns=markers, index=df.index),phenotypes

    def isolateCluster(self, col_name):
        """ Parse out cluster name from column name"""
        matches = re.findall('#\d+', col_name)
        if matches:
            return matches[0][1:]  # exclude first character!
        else:
            print('Error. Check the ce column names')
            return None




    def addTsne(self, data):
        dir = os.getcwd()
        for file in os.listdir(dir):
            if file.find('info_hal') > -1:
                for file in os.listdir(str(dir) + '/info_hal'):
                    if file.find('tsne') > -1:
                        tsne = pd.read_pickle(str(dir) + '/info_hal/' + file)
                        break
        col1 = []
        col2 = []
        for x in tsne:
            col1.append(x[0])
            col2.append(x[1])
        data["TSNE1"] = col1
        data["TSNE2"] = col2

        return data
    
     
        
    def assignRegion(self,number,dist):
 
        for i,value in enumerate(self.phenotype_complexity):
            if number < np.percentile(dist,value*100):
                return i+1
        return len(self.phenotype_complexity)


    def add_phenotypes(self,data):
        phens = pd.DataFrame()
        col_nums = []
        for i, col in enumerate(data.columns):  # iterate through clusters
            try:
                num = int(col)
            except ValueError:
                num = int(self.isolateCluster(col))
            col_nums.append(num)
        phenotype_dict = {col_num: i for i, col_num in enumerate(set(col_nums))}
        for i, num in enumerate(col_nums):
            phenotype = self.phenotypes[phenotype_dict[num]]
            phens[data.columns[i]] = [phenotype]
        phens.index = ["Phenotypes"]
        data = data.append(phens)
        return data

