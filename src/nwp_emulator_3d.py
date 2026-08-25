import torch
import torch.nn as nn

class NWPEmulator3D(nn.Module):
    def __init__(self):
        super(NWPEmulator3D, self).__init__()
        
        # Input: 5 channels (Temp, RH, U, V, HGT), 15 vertical levels
        self.encoder = nn.Sequential(
            nn.Conv3d(in_channels=5, out_channels=32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv3d(in_channels=32, out_channels=64, kernel_size=3, padding=1),
            nn.ReLU()
        )
        
        # Decoder: Maps back to the 5 atmospheric variables for Time T+1
        self.decoder = nn.Sequential(
            nn.Conv3d(in_channels=64, out_channels=32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv3d(in_channels=32, out_channels=5, kernel_size=3, padding=1)
        )

    def forward(self, x):
        features = self.encoder(x)
        out = self.decoder(features) 
        return out