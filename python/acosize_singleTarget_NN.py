import torch
from torch import nn
from torch.utils.data import DataLoader
from sklearn.preprocessing import StandardScaler
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

species = 'SAI'

# --------------------------------------
# define classes and functions
# --------------------------------------
class AcoSizeDataset(torch.utils.data.Dataset):
    def __init__(self, x, y, scale_data=True):
        if not torch.is_tensor(x) and not torch.is_tensor(y):
            # Apply scaling if necessary
            if scale_data:
                x = StandardScaler().fit_transform(x)
            self.x = torch.from_numpy(x)
            self.y = torch.from_numpy(y)

    def __len__(self):
        return len(self.x)

    def __getitem__(self, i):
        return self.x[i], self.y[i]


# Fully connected neural network with one hidden layer
class NeuralNet(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        super(NeuralNet, self).__init__()
        self.l1 = nn.Linear(input_size, hidden_size)  # first layer
        self.relu = nn.ReLU()  # activation function
        self.l2 = nn.Linear(hidden_size, num_classes)  # second layer

    def forward(self, x):
        out = self.l1(x)
        out = self.relu(out)
        out = self.l2(out)
        # no activation and no softmax at the end
        return out

# --------------------------------------
# get things going
# --------------------------------------
raw_df = pd.read_csv('G:/acosize_ts/results/acosize_200_'+species+'_features and targets.csv',
                     sep=",")

# data = raw_df[['angle',
#                'SNR',
#                'meanVal',
#                'sdVal',
#                'PC1',
#                'PC2']].values
data = raw_df[['angle',
               'SNR',
               'meanVal',
               'sdVal',
               'PC1',
               'PC2',
               'nullsN',
               'nullsProm',
               'mu',
               'sigma',
               'skew',
               'kurt',
               'Npeaks',
               'promPeaks']].values

nSamples = np.shape(data)[0]
input_size = np.shape(data)[1]
hidden_size = 150 # number of nodes in hidden layer
num_classes = np.shape(np.unique(raw_df[['lengthGroupNum']].values))[0] # number of classes, 0-9
num_epochs = 50 # number of times we go through the entire dataset
batch_size = 50 # number of samples in one forward/backward pass
learning_rate = 1e-4 # learning rate

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = NeuralNet(input_size, hidden_size, num_classes).to(device)

target = raw_df[['lengthGroupNum']].values
mapper, ind = np.unique(target, return_inverse=True)
target = ind

# sample_idx = torch.randint(len(data), size=(nSamples,1))
# sample_idx = np.reshape(sample_idx,[-1])
# sampleTrain = sample_idx[1:trainingSize]
# sampleTest = sample_idx[trainingSize+1:trainingSize+testSize+1]
#
# dataTraining = AcoSizeDataset(data[sampleTrain,], target[sampleTrain])
# dataTest = AcoSizeDataset(data[sampleTest,], target[sampleTest])

train_loader = torch.utils.data.DataLoader(dataset=AcoSizeDataset(data, target), batch_size=batch_size,
                                           shuffle=True)

test_loader = torch.utils.data.DataLoader(dataset=AcoSizeDataset(data, target), batch_size=batch_size,
                                          shuffle=False)

# loss and optimizer
criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(),lr=learning_rate)

# training loop
n_total_steps = len(train_loader) # number of batches in training set

for epoch in range(num_epochs):
    for i, data in enumerate(train_loader):

        # Get and prepare inputs
        inputs, targets = data
        inputs, targets = inputs.float(), targets.long()
        #targets = targets.reshape((targets.shape[0], 1))

        # forward pass
        outputs = model(inputs)
        loss = criterion(outputs,targets)

        # backward pass
        optimizer.zero_grad() # set gradients to zero
        loss.backward()     # backpropagation
        optimizer.step()   # update weights

        if (i+1) % 100 == 0: # print every 100 steps
            print(f'epoch {epoch+1}/{num_epochs}, step {i+1}/{n_total_steps}, loss = {loss.item():.4f}')

print('Finished training')

# test the model
with torch.no_grad(): # we don't need gradients in the testing phase
    n_correct = 0
    n_samples = 0
    for inputs, targets in test_loader:

        # Get and prepare inputs
        inputs, targets = inputs.float(), targets.long()

        outputs = model(inputs)                     # 100,10

        # value, index
        _, predictions = torch.max(outputs,1) # 1 is the dimension
        n_samples += targets.shape[0] # number of samples in the current batch
        n_correct += (predictions == targets).sum().item()  # number of correct predictions

    acc = 100.0 * n_correct / n_samples  # accuracy
    print(f'accuracy = {acc}')


