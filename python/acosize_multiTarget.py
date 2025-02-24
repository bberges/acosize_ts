import torch
from torch import nn
from torch.utils.data import DataLoader
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import pi


# --------------------------------------
# define classes and functions
# --------------------------------------
class AcoSizeDataset(torch.utils.data.Dataset):
    """
  Prepare the Boston dataset for regression
  """

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


class MLP(nn.Module):
    """
    Multilayer Perceptron for regression.
  """

    def __init__(self):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Linear(5, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 2)
        )

    def forward(self, x):
        """
      Forward pass
    """
        return self.layers(x)


def accuracy(model, ds, pct_close):
    n_correct_L = 0
    n_wrong_L = 0
    n_correct_A = 0
    n_wrong_A = 0
    for i in range(len(ds)):
        x = ds[i][0].float()  # 2-d inputs
        y = ds[i][1].float()  # 2-d target
        with torch.no_grad():
            oupt = model(x)  # computed income

        if torch.abs(oupt[0] - y[0]) < torch.abs(pct_close * y[0]):
            n_correct_L += 1
        else:
            n_wrong_L += 1

        if torch.abs(oupt[1] - y[1]) < torch.abs(pct_close * y[1]):
            n_correct_A += 1
        else:
            n_wrong_A += 1
    acc_L = (n_correct_L * 1.0) / (n_correct_L + n_wrong_L)
    acc_A = (n_correct_A * 1.0) / (n_correct_A + n_wrong_A)
    return acc_L,acc_A


def angle(angle1, angle2):
    a = angle1 - angle2
    a = (a + 180) % 360 - 180
    print(abs(a))


# --------------------------------------
# get things going
# --------------------------------------
raw_df = pd.read_csv('G:/acosize_ts/results/acosize_200_SAI_features and targets.csv',
                     sep=",")


# plt.hist(raw_df.loc[:,'angle'])
# plt.show()
#
# plt.plot(raw_df.loc[:,'angle'])
# plt.show()
#
# plt.plot(np.log10(np.abs(np.tan(raw_df.loc[:,'angle']*2*pi/360))))
# plt.show()
#
# plt.scatter(raw_df.loc[:,'angle'],np.log10(np.abs(np.tan(raw_df.loc[:,'angle']*2*pi/360))))
# plt.show()

raw_df.loc[:,'angle'] = np.abs(np.tan(raw_df.loc[:,'angle']*2*pi/360))
raw_df.loc[:,'angle'] = np.arctan(raw_df.loc[:,'angle'])
#raw_df.loc[:,'angle'] = np.abs(np.tan(raw_df.loc[:,'angle']))
#raw_df.loc[raw_df['angle'] > 180,'angle'] = raw_df.loc[raw_df['angle'] > 180,'angle'] -180
#
# plt.hist(raw_df.loc[:,'angle'])
# plt.show()

data = raw_df[['SNR', 'meanVal', 'sdVal', 'PC1', 'PC2']].values
target = raw_df[['length','angle']].values

dataset = AcoSizeDataset(data, target)

trainloader = torch.utils.data.DataLoader(dataset, batch_size=100, shuffle=True)  # , shuffle=True, num_workers=1)

# Initialize the MLP
mlp = MLP()

# Define the loss function and optimizer
loss_function = nn.L1Loss()
optimizer = torch.optim.Adam(mlp.parameters(), lr=1e-4)

losses = []
# Run the training loop
for epoch in range(0, 20):  # 5 epochs at maximum

    # Print epoch
    print(f'Starting epoch {epoch + 1}')

    # Set current loss value
    current_loss = 0.0

    # Iterate over the DataLoader for training data
    for i, data in enumerate(trainloader, 0):

        # Get and prepare inputs
        inputs, targets = data
        inputs, targets = inputs.float(), targets.float()
        targets = targets.reshape((targets.shape[0], 2))

        # Zero the gradients
        optimizer.zero_grad()

        # Perform forward pass
        outputs = mlp(inputs)

        # Compute loss
        loss = loss_function(outputs, targets)

        # Perform backward pass
        loss.backward()

        # Perform optimization
        optimizer.step()

        losses.append(loss.item())

        # Print statistics
        current_loss += loss.item()
        if i % 10 == 0:
            print('Loss after mini-batch %5d: %.3f' %
                  (i + 1, current_loss / 500))
            current_loss = 0.0

# Process is complete.
print('Training process has finished.')

# plt.plot(losses)
# plt.xlabel("no. of iterations")
# plt.ylabel("total loss")
# plt.show()

minAngle = 0
maxAngle = 180
filt_df = raw_df[(raw_df['angle'] <= maxAngle) & (raw_df['angle'] >= minAngle)]
data_test = filt_df[['SNR', 'meanVal', 'sdVal', 'PC1', 'PC2']].values
target_test = filt_df[['length','angle']].values
np.shape(target_test)

dataset_test = AcoSizeDataset(data_test, target_test)

print(accuracy(mlp, dataset_test, 0.1))

oupt = np.empty((2,len(dataset)))
for i in range(len(dataset)):
    x = dataset[i][0].float()  # 2-d inputs
    y = dataset[i][1].float()  # 2-d target
    with torch.no_grad():
        oupt[:,i] = mlp(x).detach().numpy() # prediction

print(accuracy(mlp, dataset, 0.3))

plt.scatter(np.abs(oupt[1,:])*360/2/pi, target[:,1]*360/2/pi)
plt.show()
plt.close()

plt.scatter(oupt[0,:], target[:,0])
plt.show()
plt.close()

print(y)
print(oupt)

angle(y[1],oupt[1])