import os
import sys
import random
import json
import numpy as np
import torch
import argparse
from model_dnn import *
import matplotlib.pyplot as plt
import sklearn.metrics
import sklearn.model_selection

TIRAMISU_DEBUG = False

def get_params(params_fpath):
    
    if params_fpath is None:
        params = {
            'dropout': 0.5,
            'n_hidden': [21, 40],
            'activate_func': 'relu'
        }
    
    else:
        with open(params_fpath, 'r') as jsonfh:
            params = json.load(jsonfh)
        
        n_hidden = []
        if 'n_hidden_0' in params:
            n_hidden.append(params['n_hidden_0'])
            del params['n_hidden_0']
        if 'n_hidden_1' in params:
            n_hidden.append(params['n_hidden_1'])
            del params['n_hidden_1']
        if 'n_hidden_2' in params:
            n_hidden.append(params['n_hidden_2'])
            del params['n_hidden_2']
    
        params['n_hidden'] = n_hidden
    
    return params
    



def train_cv(model, weight, feature_type, cv_dataset, batch_size=1024, epochs=50, n_tries=5):
    
    if weight is None:
        raise ValueError('argument `weight` cannot be None.')
    
    # set up candidates of hyper-parameters of network architecture
    if model == 'L1':
        dropouts = [0, 0.5]
        n_hiddens = [[i] for i in [8, 12, 16, 20, 24, 28, 32, 36, 40]]
        activate_funcs = ['relu']
        if TIRAMISU_DEBUG:
            dropouts = [0, 0.5]
            n_hiddens = [[i] for i in [8, 22]]
            activate_funcs = ['relu']
        
    elif model == 'L2':
        dropouts = [0.5]
        n_hiddens = [[i, j] for i in [8, 16, 24, 32, 40]
                            for j in [8, 16, 24, 32, 40]]
        activate_funcs = ['relu']
        if TIRAMISU_DEBUG:
            dropouts = [0, 0.5]
            n_hiddens = [[i, j] for i in [ 8, 22]
                                for j in [18, 24]]
            activate_funcs = ['relu']
        
    elif model == 'L3':
        dropouts = [0.5]
        n_hiddens = [[i, j, k] for i in [ 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
                               for j in [16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40]
                               for k in [20, 24, 28, 32, 36, 40]]
        activate_funcs = ['relu']
    
    
    # 10-fold cv to evaluate model
    val_outputs = None
    for i in range(1, 11):
        train_dataset = pd.read_csv(os.path.join(cv_dataset, 'train_cv.' + str(i) + '.std.csv'), header=0, index_col=0)
        valid_dataset = pd.read_csv(os.path.join(cv_dataset, 'valid_cv.' + str(i) + '.std.csv'), header=0, index_col=0)
        
        # glid search / 2-fold cv to determine hyper-params of DNN
        seed = abs(hash(cv_dataset)) % (10 ** 8)
        best_loss = 1e100 * 1.0
        best_params = None
        best_weight = None
        for activate_func in activate_funcs:
            for n_hidden in n_hiddens:
                for dropout in dropouts:
                    random.seed(seed)
                    np.random.seed(seed)
                    torch.manual_seed(seed)
 
                    weight_ = '{}_cv.{}_{}_{}_{}.pth'.format(weight[:-4], i, activate_func, '-'.join([str(i) for i in n_hidden]), dropout)
                    param_ = {'dropout': dropout, 'n_hidden': n_hidden, 'activate_func': activate_func}
                    loss_cv_ = []
                    kf = sklearn.model_selection.KFold(n_splits=2, random_state=seed, shuffle=True)
                    print(param_)
                    for train_index, valid_index in kf.split(train_dataset):
                        _train = train_dataset.iloc[train_index, :]
                        _valid = train_dataset.iloc[valid_index, :]
                        # train
                        loss_hist = train(weight_, feature_type, _train, _valid, batch_size, epochs, param_, False, n_tries)
                        loss_cv_.append(sum(loss_hist['valid']) / len(loss_hist['valid']))  # average of `n_tries`
                    loss_cv_mean_ = sum(loss_cv_) / len(loss_cv_)
                    if loss_cv_mean_ < best_loss:
                        best_loss = loss_cv_mean_
                        best_params = param_
                        best_weight = weight_
        
        loss_hist = train(best_weight, feature_type, train_dataset, valid_dataset, batch_size, epochs, best_params, True, n_tries)
        pred_values = validate(best_weight, feature_type, best_params, valid_dataset, output=None)
        val_outputs_ = pd.DataFrame({'true': valid_dataset.loc[:, 'incidence'].values, 'cv': i, 'predicted': pred_values})

        val_outputs = pd.concat([val_outputs, val_outputs_], axis=0)
    
    
    print(val_outputs)
    val_outputs.to_csv(weight[:-4] + '_cvresults.tsv', header=True, index=False, sep='\t')





def train(weight, feature_type, train_dataset, valid_dataset, batch_size, epochs, params = None, save_weight=True, n_try=10):
    if TIRAMISU_DEBUG:
        epochs = 3
        n_try = 2
    
    if params is None:
        params = get_params()
    
    if isinstance(params['n_hidden'], int):
            params['n_hidden'] = [params['n_hidden']]
    
    if weight is None:
        raise ValueError('argument `weight` cannot be None.')
    
    best_loss_ = 1e10 * 1.0
    best_model_ = None
    min_loss = {'train': [], 'valid': []}
    train_history_ = None
    
    # since little dataset causes unstable models,
    # try `n_try` times to create model and get the best one according to validation scores
    for i in range(n_try):
        jppnet = Jppnet(feature_type=feature_type, n_hidden=params['n_hidden'], dropout=params['dropout'], activate_func=params['activate_func'])
        train_history = jppnet.train(train_dataset, valid_dataset, batch_size, epochs)
        print('try {} -- train loss: {};  valid loss: {}'.format(
                i, train_history.loc[:, 'train'].values[-1],
                   train_history.loc[:, 'valid'].values[-1]))

        min_loss['train'].append(train_history.loc[:, 'train'].values[-1])
        min_loss['valid'].append(train_history.loc[:, 'valid'].values[-1])
        if train_history.loc[:, 'valid'].values[-1] < best_loss_:
            best_loss_  = train_history.loc[:, 'valid'].values[-1]
            best_model_ = jppnet
            train_history_ = train_history
    
    print('mean  -- train loss mean: {};  valid loss mean: {}'.format(
           np.mean(min_loss['train']), np.mean(min_loss['valid'])))
    
    if save_weight:
        jppnet.save(weight)
        # loss of `n_try` tries
        pd.DataFrame(min_loss).to_csv(weight + '_stats.tsv', index=False, header=True, sep='\t')
        # train history of the best model
        pd.DataFrame(train_history_).to_csv(weight + '_trainhistory.tsv', index=False, header=True, sep='\t')
    
    return min_loss




def validate(weight, feature_type, params, valid_dataset, output=None):
    
    jppnet = Jppnet(feature_type=feature_type, n_hidden=params['n_hidden'], dropout=params['dropout'], activate_func=params['activate_func'])
    jppnet.load_weights(weight)
    predicted_df = jppnet.validate(valid_dataset)

    if output is not None:
        predicted_df.to_csv(output + '.tsv', index=False, header=True, sep='\t')
        
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.scatter(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted'])
        ax.set_xlabel('label')
        ax.set_ylabel('predicted')
        fig.savefig(output + '.png')
    
    print('--- validation results ---')
    print('R2 :   {}'.format(sklearn.metrics.r2_score(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted'])))
    print('MSE:   {}'.format(sklearn.metrics.mean_squared_error(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted'])))
    print('RMSE:  {}'.format(np.sqrt(sklearn.metrics.mean_squared_error(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted']))))
    print('MAE:   {}'.format(sklearn.metrics.mean_absolute_error(predicted_df.loc[:, 'label'], predicted_df.loc[:, 'predicted'])))
    print('--------------------------')
    
    return predicted_df.loc[:, 'predicted'].values
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Let dragonfly fly!')
     
    parser.add_argument('--model', default='L1')
    parser.add_argument('--feature-type', default='category')
    parser.add_argument('--weight', default='./weights/test.pth')
    parser.add_argument('--output', default='./weights/test.txt')
    parser.add_argument('--cv-dataset', default=None)
    parser.add_argument('--train-dataset', default='../formatted_data/test/train_std.tsv')
    parser.add_argument('--valid-dataset', default='../formatted_data/test/valid_std.tsv')
    parser.add_argument('--epochs', default=50, type=int)
    parser.add_argument('--batch-size', default=1024, type=int)
    parser.add_argument('--mode', default='train')
    parser.add_argument('--params', default=None)
    
    args = parser.parse_args()
    
    if args.mode == 'cv':
        train_cv(args.model, args.weight, args.feature_type, args.cv_dataset, args.batch_size, args.epochs)
        
    elif args.mode == 'train':
        params = get_params(args.params)
        train(args.weight, args.feature_type, args.train_dataset, args.valid_dataset, args.batch_size, args.epochs, params, True, 5)
    
    elif args.mode == 'valid':
        params = get_params(args.params)
        validate(args.weight, args.feature_type, params, args.valid_dataset, args.output)
    
    elif args.mode == 'inference':
        pass    
    
    
    


