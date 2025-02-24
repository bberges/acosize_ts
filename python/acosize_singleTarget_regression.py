import torch
from torch import nn
from torch.utils.data import DataLoader
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

species = 'MAC'

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

    def __init__(self,input_size):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Linear(input_size, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 1)
        )

    def forward(self, x):
        """
      Forward pass
    """
        return self.layers(x)


def accuracy(model, ds, pct_close):
    n_correct = 0
    n_wrong = 0
    for i in range(len(ds)):
        x = ds[i][0].float()  # 2-d inputs
        y = ds[i][1].float()  # 2-d target
        with torch.no_grad():
            oupt = model(x)  # computed income

        if torch.abs(oupt - y) < torch.abs(pct_close * y):
            n_correct += 1
        else:
            n_wrong += 1
    acc = (n_correct * 1.0) / (n_correct + n_wrong)
    return acc


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

target = raw_df[['length']].values

dataset = AcoSizeDataset(data, target)

trainloader = torch.utils.data.DataLoader(dataset, batch_size=50, shuffle=True)  # , shuffle=True, num_workers=1)

# Initialize the MLP
mlp = MLP(np.shape(data)[1])

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
        targets = targets.reshape((targets.shape[0], 1))

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

plt.plot(losses)
plt.xlabel("no. of iterations")
plt.ylabel("total loss")
plt.savefig('loss_profile'+species+'.png', dpi=300)
plt.close()


oupt = np.empty((1,len(dataset)))
for i in range(len(dataset)):
    x = dataset[i][0].float()  # 2-d inputs
    y = dataset[i][1].float()  # 2-d target
    with torch.no_grad():
        oupt[:,i] = mlp(x).detach().numpy() # prediction

outFrame =np.column_stack((raw_df[['angle']].to_numpy(),target,np.transpose(oupt)))
pd.DataFrame(data=outFrame,columns=['angle','length','length_pred']).to_csv('G:/acosize_ts/results/pred_'+species+'.csv', index=False)
#
# plt.scatter(oupt[0,:], target[:,0])
# plt.ylabel('targets (cm)')
# plt.xlabel('predictions (cm)')
# plt.title('SAI')
# plt.savefig('pred_SAI.png', dpi=300)
# plt.close()
#
# minAngle = 240
# maxAngle = 280
# filt_df = raw_df[(raw_df['angle'] <= maxAngle) & (raw_df['angle'] >= minAngle)]
# data_test = filt_df[['angle', 'SNR', 'meanVal', 'sdVal', 'PC1', 'PC2']].values
# target_test = filt_df[['length']].values
# np.shape(target_test)
#
# dataset_test = AcoSizeDataset(data_test, target_test)
#
# accuracy(mlp, dataset_test, 0.3)
