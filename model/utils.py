import os
import sys
import random
import numpy as np
import pandas as pd



PREF2ID = {'hokkaido' :  0, 'aomori'   :  1, 'iwate'   :  2, 'miyagi'   :  3, 'akita'    :  4,
           'yamagata' :  5, 'fukushima':  6, 'ibaraki' :  7, 'tochigi'  :  8, 'gunma'    :  9,
           'saitama'  : 10, 'chiba'    : 11, 'tokyo'   : 12, 'kanagawa' : 13, 'niigata'  : 14,
           'toyama'   : 15, 'ishikawa' : 16, 'fukui'   : 17, 'yamanashi': 18, 'nagano'   : 19,
           'gifu'     : 20, 'shizuoka' : 21, 'aichi'   : 22, 'mie'      : 23, 'shiga'    : 24,
           'kyoto'    : 25, 'osaka'    : 26, 'hyogo'   : 27, 'nara'     : 28, 'wakayama' : 29,
           'tottori'  : 30, 'shimane'  : 31, 'okayama' : 32, 'hiroshima': 33, 'yamaguchi': 34,
           'tokushima': 35, 'kagawa'   : 36, 'ehime'   : 37, 'kochi'    : 38, 'fukuoka'  : 39,
           'saga'     : 40, 'nagasaki' : 41, 'kumamoto': 42, 'oita'     : 43, 'miyazaki' : 44,
           'kagoshima': 45, 'okinawa'  : 46}




def load_dataset(fpath):
    df = pd.read_csv(fpath, header=0, index_col=0)
    return df




def randomize_dataset(df, method='shuffle', seed=0):
    np.random.seed(seed)
    df_randomized = df.copy()
    
    if method == 'shuffle':
        df_randomized = df_randomized.sample(frac=1, replace=False, random_state=seed)
    elif method == 'randomize':
        for _ in range(df_randomized.shape[1]):
            np.random.seed(seed + _ * 2)
            random.seed(seed + _ * 2)
            df_randomized.iloc[:, _] = random.sample(df_randomized.iloc[:, _].values.tolist(), len(df_randomized.iloc[:, _]))
    else:
        raise ValueError('not supported!')
    
    return df_randomized



def month2onehot(x):
    y = np.identity(12)[x.values - 1]
    y = y[:, 1:]   # drop the first column
    y = y.astype(np.int16)
    return y
    


def pref2onehot(x):
    y = x.replace(PREF2ID).astype('int64')
    y = np.identity(47)[y.values]
    y = y[:, 1:]   # drop the first column
    y = y.astype(np.int16)
    return y
    


def get_xy(df, feature_type='category'):
    subset_id = df.iloc[:, -2].values
    y = df.iloc[:, -3].values
    X = None
    
    if feature_type == 'category':
        m = month2onehot(df.loc[:, 'Month'])
        p = pref2onehot(df.loc[:, 'Pref'])
        X = np.concatenate([m, p], 1)
    elif feature_type == 'decimal':
        X = df.loc[:, ['longi', 'lati', 'airtemp', 'precip']].dropna(how='any', axis=0).values
    else:
        raise ValueError('set `category` or `decimal` in `get_xy` function.')
    
    
    return X, y, subset_id







