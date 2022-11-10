import os
import numpy as np
import pandas as pd
import torch
import torch.optim
import torch.nn
import torch.nn.functional as F
import torchvision
import torch.autograd

# os.environ['KMP_DUPLICATE_LIB_OK']='True'


PREFECTURE_DICT = {
                   'hokkaido': 1,
                   'aomori' : 2,
                   'iwate' : 3,
                   'miyagi' : 4,
                   'akita' : 5,
                   'yamagata' : 6,
                   'fukushima': 7,
                   'ibaraki': 8,
                   'tochigi': 9,
                   'gunma': 10,
                   'saitama': 11,
                   'chiba': 12,
                   'tokyo': 13,
                   'kanagawa': 14,
                   'niigata': 15,
                   'toyama': 16,
                   'ishikawa': 17,
                   'fukui': 18,
                   'yamanashi': 19,
                   'nagano': 20,
                   'gifu': 21,
                   'shizuoka': 22,
                   'aichi': 23,
                   'mie': 24,
                   'shiga': 25,
                   'kyoto': 26,
                   'osaka': 27,
                   'hyogo': 28,
                   'nara': 29,
                   'wakayama': 30,
                   'tottori': 31,
                   'shimane': 32,
                   'okayama': 33,
                   'hiroshima': 34,
                   'yamaguchi': 35,
                   'tokushima': 36,
                   'kagawa': 37,
                   'ehime': 38,
                   'kochi': 39,
                   'fukuoka': 40,
                   'saga': 41,
                   'nagasaki': 42,
                   'kumamoto': 43, 
                   'oita': 44,
                   'miyazaki': 45,
                   'kagoshima': 46,
                   'okinawa': 47
                  }







class JppnetArch(torch.nn.Module):
    
    def __init__(self, n_input=None, n_hidden=12, dropout=0, activate_func='relu'):
        super(JppnetArch, self).__init__()
        
        if isinstance(n_hidden, int):
            n_hidden = [n_hidden]
        
        # input layer
        self.input = torch.nn.Linear(n_input, n_hidden[0])

        # hidden layer
        hidden_layers = []
        if len(n_hidden) > 1:
            for i in range(1, len(n_hidden)):
                if activate_func == 'relu':
                    hidden_layers += [torch.nn.ReLU(inplace=True)]
                else:
                    hidden_layers += [torch.nn.Sigmoid()]
                if dropout > 0:
                    hidden_layers += [torch.nn.Dropout(p=dropout)]
                hidden_layers += [torch.nn.Linear(n_hidden[i - 1], n_hidden[i])]

        if activate_func == 'relu':
            hidden_layers += [torch.nn.ReLU(inplace=True)]
        else:
            hidden_layers += [torch.nn.Sigmoid()]
        if dropout > 0:
            hidden_layers += [torch.nn.Dropout(p=dropout)]
        self.fc = torch.nn.Sequential(*hidden_layers)

        # output layer
        self.output = torch.nn.Linear(n_hidden[-1], 1)
    
    
    def forward(self, x):
        x = self.input(x)
        x = self.fc(x)
        x = self.output(x)
        return x






class textDataset(torch.utils.data.Dataset):
    
    
    def __init__(self, dat_fpath, feature_type = 'category', is_test=False):
        self.is_test = is_test
        
        if (type(dat_fpath) is str):
            self.X, self.y = self.__read_csv(dat_fpath, feature_type)
        else:
            # pandas data.frame
            self.X, self.y =  self.__read_df(dat_fpath, feature_type)
        
        self.size = len(self.y)
    
    def __read_csv(self, dat_fpath, feature_type):
        '''
        load dataset from CSV file
        '''
        
        X = []
        y = []
        
        with open(dat_fpath, 'r') as datfh:
            for dat_buf in datfh:
                _i, year, month, pref, season, incidence, lng, lat, temp, rain, area = dat_buf.replace('\n', '').replace('"', '').split(',')
                
                if _i == '':
                    continue  # data header
                
                _x = []
                _y = incidence
                
                # features
                if feature_type == 'category':
                    # use one-hot encoding for month and prefecture
                    _x.extend(self.__month2onehot(month))
                    _x.extend(self.__prefecture2onehot(pref))
                elif feature_type == 'decimal':
                    # use temperature and precipication instead of month and prefecture name
                    _x = [float(temp), float(rain), float(lat), float(lng)]
                else:
                    raise ValueError('set `category` or `decimal` in `get_xy` function.')
                
                # labels
                if self.is_test:
                    X.append(_x)
                    y.append('NA')
                else:
                    if _y != 'NA':
                        _y = [float(_y)]
                        #_y = [np.log10(float(_y) + 1)]
                        X.append(_x)
                        y.append(_y)
        
        X = torch.from_numpy(np.array(X)).float()
        y = torch.from_numpy(np.array(y)).float()
        
        return X, y
        
    
    def __read_df(self, df, feature_type):
        '''
        arrange dataset from data.frame
        '''
        
        X = []
        y = []
        
        year = df.loc[:, 'Year'].values
        month = df.loc[:, 'Month'].values
        pref = df.loc[:, 'Pref'].values
        incidence = df.loc[:, 'incidence'].values
        lng = df.loc[:, 'longi'].values
        lat = df.loc[:, 'lati'].values
        temp = df.loc[:, 'airtemp'].values
        rain = df.loc[:, 'precip'].values
        
        for i in range(len(month)):
            _x = []
            _y = incidence[i]
                
            # features
            if feature_type == 'category':
                # use one-hot encoding for month and prefecture
                _x.extend(self.__month2onehot(month[i]))
                _x.extend(self.__prefecture2onehot(pref[i]))
            elif feature_type == 'decimal':
                # use temperature and precipication instead of month and prefecture name
                _x = [float(temp[i]), float(rain[i]), float(lat[i]), float(lng[i])]
            else:
                raise ValueError('set `category` or `decimal` in `get_xy` function.')
                
            # labels
            if self.is_test:
                X.append(_x)
                y.append('NA')
            else:
                if _y != 'NA':
                    _y = [float(_y)]
                    X.append(_x)
                    y.append(_y)
        
        X = torch.from_numpy(np.array(X)).float()
        y = torch.from_numpy(np.array(y)).float()
        
        return X, y
        
 

    def __len__(self):
        return self.y.shape[0]
    
    
    def __getitem__(self, i):
        _x = self.X[i]
        _y = self.y[i]
        
        if self.is_test:
            return _x
        else:
            return _x, _y
    
    
    def __month2onehot(self, m):
        m_vec = np.zeros(12)
        m_vec[int(m) - 1] = 1
        return m_vec.tolist()
        #m = int(m)
        #if m == 12 or m == 1 or m == 2:
        #    m_vec = [1, 0, 0, 0]
        #elif m == 3 or m == 4 or m == 5:
        #    m_vec = [0, 1, 0, 0]
        #elif m == 6 or m == 7 or m == 8:
        #    m_vec = [0, 0, 1, 0]
        #elif m == 9 or m == 10 or m == 11:
        #    m_vec = [0, 0, 0, 1]
        #elif m == 9 or m == 10 or m == 11:
        #    m_vec = [0, 0, 0, 1]
        #elif m == 9 or m == 10 or m == 11:
        #    m_vec = [0, 0, 0, 1]
        #
        #return m_vec
    
    def __prefecture2onehot(self, p):
        p_vec = np.zeros(len(PREFECTURE_DICT))
        p_vec[PREFECTURE_DICT[p.lower()] - 1] = 1
        return p_vec.tolist()







class Jppnet():
    
    def __init__(self, feature_type=None, n_hidden=12, dropout=0.5, activate_func='relu'):
        
        if feature_type == 'category':
            n_input = 12 + 47
        elif feature_type == 'decimal':
            n_input = 4
        
        self.model = JppnetArch(n_input=n_input, n_hidden=n_hidden, dropout=dropout, activate_func=activate_func) 
        self.__best_model = JppnetArch(n_input=n_input, n_hidden=n_hidden, dropout=dropout, activate_func=activate_func) 
        #self.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        self.device = torch.device('cpu')
        self.model.to(self.device)
        self.feature_type = feature_type
        
    def load_weights(self, weights_fpath):
        self.model.load_state_dict(torch.load(weights_fpath))
    
    
    def save(self, weights_fpath):
        torch.save(self.model.state_dict(), weights_fpath)
    
    
    def train(self, train_data, valid_data, batch_size=1024, num_epochs=100, num_workers=4):
        # load train and validation dataset
        datasets = {
            'train': textDataset(train_data, self.feature_type),
            'valid': textDataset(valid_data, self.feature_type)
        }
        dataloader = {
            'train': torch.utils.data.DataLoader(datasets['train'],
                                                  batch_size=batch_size, shuffle=True, num_workers=num_workers),
            'valid': torch.utils.data.DataLoader(datasets['valid'],
                                                  batch_size=batch_size, num_workers=num_workers)
        }
        dataset_sizes = {
            'train': len(datasets['train']),
            'valid': len(datasets['valid']),
        }
        
        # training parameters
        criterion = torch.nn.MSELoss()
        optimizer = torch.optim.Adam(self.model.parameters(), lr=1e-2)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=20, gamma=0.1)
        
        # training and validation
        epoch_loss = {'train': [], 'valid': []}
        _last_epoch_loss = 1e10
        _last_epoch_loss_count = 0
        for epoch in range(num_epochs):
            
            for phase in ['train', 'valid']:

                if phase == 'train':
                    self.model.train()
                else:
                    self.model.eval()
            
                running_loss = 0.0
                running_corrects = 0
        
                for inputs, labels in dataloader[phase]:
                    inputs = inputs.to(self.device)
                    labels = labels.to(self.device)
                    
                    # forward
                    with torch.set_grad_enabled(phase == 'train'):
                        outputs = self.model(inputs)
                        loss = criterion(outputs, labels)

                        if phase == 'train':
                            optimizer.zero_grad()
                            loss.backward()
                            optimizer.step()

                    running_loss += loss.item() * inputs.size(0)

                if phase == 'train':
                    scheduler.step()

                _epoch_loss = running_loss / dataset_sizes[phase]
                
                if phase == 'valid':
                    if _epoch_loss < _last_epoch_loss:
                        _last_epoch_loss = _epoch_loss
                        _last_epoch_loss_count = 0
                        self.__best_model.load_state_dict(self.model.state_dict())
                    else:
                        _last_epoch_loss_count += 1
            
                #print('{} Loss: {:.4f}'.format(phase, _epoch_loss))
                epoch_loss[phase].append(_epoch_loss)
            
            if _last_epoch_loss_count > 20:
                break
        
        epoch_loss = pd.DataFrame(epoch_loss)
        self.model.load_state_dict(self.__best_model.state_dict())
        return epoch_loss
    
    
    
    
    def validate(self, test_data, batch_size=None):
        
        self.model.eval()
        
        dataset = textDataset(test_data, self.feature_type)
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=1024, shuffle=False)
        
        pred = {'label': [], 'predicted': []}
        
        with torch.set_grad_enabled(False):
            for inputs, labels in dataloader:
                inputs = inputs.to(self.device)
                pred['predicted'].extend(self.model(inputs).data.cpu().numpy().flatten())
                pred['label'].extend(labels.data.cpu().numpy().flatten())
        
        pred = pd.DataFrame(pred)
        
        return pred
    
    
    
    def inference(self, test_data, batch_size=None):
        
        self.model.eval()
        
        dataset = textDataset(test_data, self.feature_type, is_test)
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=1024, shuffle=False)
        
        pred = {'label': [], 'predicted': []}
        
        with torch.set_grad_enabled(False):
            for inputs, labels in dataloader:
                inputs = inputs.to(self.device)
                pred['predicted'].extend(self.model(inputs).data.cpu().numpy().flatten())
                pred['label'].extend(labels.data.cpu().numpy().flatten())
        
        pred = pd.DataFrame(pred)
        
        return pred
    
    
    
    

