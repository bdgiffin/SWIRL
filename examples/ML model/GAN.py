import os
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import pandas as pd

# Load your data
filename = 'quo_data.csv'
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")


data = pd.read_csv(filename, header=None)
columns = int((data.shape[1]-3)/2)

ID = pd.to_numeric(data.iloc[1:, 0], errors='coerce').to_numpy()
D = pd.to_numeric(data.iloc[1:, 1], errors='coerce').to_numpy()
M = pd.to_numeric(data.iloc[1:, 2], errors='coerce').to_numpy()
H = data.iloc[1:, 3:columns+3].apply(pd.to_numeric, errors='coerce').to_numpy()
R = data.iloc[1:, columns+3:columns*2+3].apply(pd.to_numeric, errors='coerce').to_numpy()
H = np.abs(H)

print("ID dtype:", ID.dtype)
print("D dtype:", D.dtype)
print("M dtype:", M.dtype)
print("R dtype:", R.dtype)
print("H dtype:", H.dtype)


# Log-transform data and prepare the input tensor
radial_flat = torch.log(torch.tensor(R.T.flatten(), dtype=torch.float32))
height_flat = torch.log(torch.tensor(H.T.flatten(), dtype=torch.float32))
mass_repeated = torch.tensor(M, dtype=torch.float32).repeat(1, R.shape[1]).reshape(-1)
diameter_repeated = torch.tensor(D, dtype=torch.float32).repeat(1, R.shape[1]).reshape(-1)

log_mass = torch.log(torch.tensor(mass_repeated, dtype=torch.float32))
log_diameter = torch.log(torch.tensor(diameter_repeated, dtype=torch.float32))

input_data = torch.stack([log_mass, log_diameter, radial_flat, height_flat], dim=1)
indices = torch.randperm(input_data.size(0))
input_data = input_data[indices]

# input_data = input_data[input_data[:, 0].argsort()]  # First sort by the first column
# input_data = input_data[input_data[:, 1].argsort()]  # Then sort by the second column

# print(np.exp(input_data[200:230,:]))
# Plotting the generated data (e.g., radius vs height)
# plt.scatter(np.exp(input_data[0:200, 2]), np.exp(input_data[0:200, 3]), c='blue', label='Generated Data')
# plt.xlabel('Radius')
# plt.ylabel('Height')
# plt.title('Generated Data (Radius vs Height)')
# for i in range(np.exp(input_data[0:200, 2]).shape[0]):
#     plt.text(np.exp(input_data[i, 2]), np.exp(input_data[i, 3]), f'M:{np.exp(input_data[i, 0]):.2f}, Dia:{np.exp(input_data[i, 1]):.2f} ', fontsize=8, color='red')
# plt.legend()
# plt.show()


# Normalize the data
scaler = StandardScaler()
input_data_normalized = scaler.fit_transform(input_data)
data_tensor = torch.tensor(input_data_normalized, dtype=torch.float32)


# GAN Parameters
latent_dim = 10
data_dim = data_tensor.shape[1]
hidden_dim = 256
batch_size = 200
epochs = 5000
lr = 0.0003

# Define Generator
class Generator(nn.Module):
    def __init__(self):
        super(Generator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, data_dim),
        )
        
    def forward(self, z):
        return self.model(z)

# Define Discriminator
class Discriminator(nn.Module):
    def __init__(self):
        super(Discriminator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(data_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, 1),
            nn.Sigmoid(),
        )
        
    def forward(self, x):
        return self.model(x)

# Initialize models
generator = Generator()
discriminator = Discriminator()

# Loss and optimizers
criterion = nn.BCELoss()
optimizer_g = optim.Adam(generator.parameters(), lr=lr)
optimizer_d = optim.Adam(discriminator.parameters(), lr=lr)

# Check for existing weights
generator_weights = "generator_weights.pth"
discriminator_weights = "discriminator_weights.pth"

if os.path.exists(generator_weights) and os.path.exists(discriminator_weights):
    # Load pre-trained weights
    generator.load_state_dict(torch.load(generator_weights))
    discriminator.load_state_dict(torch.load(discriminator_weights))
    print("Weights loaded. Skipping training.")
else:
    # Training loop
    for epoch in range(epochs):
        # Train discriminator
        real_data = data_tensor[torch.randint(0, len(data_tensor), (batch_size,))]
        fake_data = generator(torch.randn(batch_size, latent_dim))
        real_labels = torch.ones(batch_size, 1)
        fake_labels = torch.zeros(batch_size, 1)
        
        # Real loss
        real_output = discriminator(real_data)
        real_loss = criterion(real_output, real_labels)
        
        # Fake loss
        fake_output = discriminator(fake_data.detach())
        fake_loss = criterion(fake_output, fake_labels)
        
        # Total loss and backprop
        loss_d = real_loss + fake_loss
        optimizer_d.zero_grad()
        loss_d.backward()
        optimizer_d.step()
        
        # Train generator
        fake_data = generator(torch.randn(batch_size, latent_dim))
        fake_output = discriminator(fake_data)
        loss_g = criterion(fake_output, real_labels)
        
        optimizer_g.zero_grad()
        loss_g.backward()
        optimizer_g.step()
        
        # Print progress
        if epoch % 500 == 0:
            print(f"Epoch {epoch} | Loss D: {loss_d.item():.4f} | Loss G: {loss_g.item():.4f}")

    # Save the trained weights
    torch.save(generator.state_dict(), generator_weights)
    torch.save(discriminator.state_dict(), discriminator_weights)
    print("Training completed. Weights saved.")

# Generate synthetic data
# Set the random seed for reproducibility
torch.manual_seed(1)  # Set seed for CPU
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(1)  # Set seed for all GPUs (if available)

# Generate synthetic data inside a no_grad context
with torch.no_grad():
    synthetic_data = generator(torch.randn(1000, latent_dim)).numpy()  # Generate data using the generator


# Denormalize and transform back
synthetic_data_denormalized = scaler.inverse_transform(synthetic_data)
data_generated = np.exp(synthetic_data_denormalized)
# print(data_generated)

# Plot synthetic data
plt.scatter(data_generated[:, 2], data_generated[:, 3], c='blue', label='Synthetic Data')
for i in range(data_generated.shape[0]):
    plt.text(data_generated[i, 2], data_generated[i, 3], f'M:{data_generated[i, 0]:.2f}, Dia:{data_generated[i, 1]:.2f} ', fontsize=8, color='red')
plt.xlabel('Radius')
plt.ylabel('Height')
plt.title('Synthetic Data (GAN)')
plt.legend()
plt.show()
