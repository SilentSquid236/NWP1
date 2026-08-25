"""
Train the 3D NWP emulator to step the atmospheric state forward one hour.

Usage (from the project root):
    python src/train_autoregressive.py                 # train with config defaults
    python src/train_autoregressive.py --epochs 5      # short run
    python src/train_autoregressive.py --resume        # continue from last checkpoint
    python src/train_autoregressive.py --dry-run       # verify wiring, no training
"""

import argparse
import sys
import time
from pathlib import Path

# Allow running as `python src/train_autoregressive.py` from the project root.
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# --- Resource cap MUST be applied before torch is imported ------------------
# OpenMP/MKL read their thread limits at import time; setting them later has
# no effect and torch would happily grab every core on this shared box.
import resources

_FRACTION = None
for i, a in enumerate(sys.argv):
    if a == "--resource-fraction" and i + 1 < len(sys.argv):
        _FRACTION = sys.argv[i + 1]
    elif a.startswith("--resource-fraction="):
        _FRACTION = a.split("=", 1)[1]

RESOURCE_PLAN = resources.apply(_FRACTION)
# ----------------------------------------------------------------------------

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split

import config
from autoregressive_dataset import AutoregressiveDataset3D
from nwp_emulator_3d import NWPEmulator3D


def control_variable_loss(prediction, target, weights):
    """
    Weighted MSE that anchors to ground truth by treating measured
    temperature as the strict control variable.

    Per-channel MSE is averaged over batch, level, and grid dims, then
    scaled by the channel weights so temperature error dominates.
    """
    mse = nn.functional.mse_loss(prediction, target, reduction="none")
    channel_loss = mse.mean(dim=[0, 2, 3, 4])          # -> [C]
    return (channel_loss * weights).mean(), channel_loss.detach()


def run_epoch(model, loader, optimizer, weights, train=True):
    model.train() if train else model.eval()
    total_loss = 0.0
    per_channel = torch.zeros(config.N_CHANNELS)
    n_batches = 0

    context = torch.enable_grad() if train else torch.no_grad()
    with context:
        for state_t, state_t1 in loader:
            if train:
                optimizer.zero_grad()

            prediction = model(state_t)
            loss, ch_loss = control_variable_loss(prediction, state_t1, weights)

            if train:
                loss.backward()
                optimizer.step()

            total_loss += loss.item()
            per_channel += ch_loss.cpu()
            n_batches += 1

    if n_batches == 0:
        return float("nan"), per_channel

    return total_loss / n_batches, per_channel / n_batches


def save_checkpoint(path, model, optimizer, epoch, best_loss):
    torch.save({
        "epoch": epoch,
        "model_state": model.state_dict(),
        "optimizer_state": optimizer.state_dict(),
        "best_loss": best_loss,
        "channels": config.CHANNELS,
        "n_levels": config.N_LEVELS,
    }, path)


def main():
    p = argparse.ArgumentParser(description="Train the 3D NWP emulator.")
    p.add_argument("--epochs", type=int, default=config.EPOCHS)
    p.add_argument("--batch-size", type=int, default=config.BATCH_SIZE)
    p.add_argument("--lr", type=float, default=config.LEARNING_RATE)
    p.add_argument("--val-split", type=float, default=0.2,
                   help="Fraction of transitions held out for validation.")
    p.add_argument("--resume", action="store_true",
                   help="Resume from the latest checkpoint if one exists.")
    p.add_argument("--dry-run", action="store_true",
                   help="Build everything and run one batch, then exit.")
    p.add_argument("--resource-fraction", type=float, default=None,
                   help="Fraction of server cores to use. Default 0.5 (half). "
                        "Parsed before torch import; shown here for --help.")
    args = p.parse_args()

    # Enforce the cap inside torch as well as via the OMP env vars.
    torch.set_num_threads(RESOURCE_PLAN["torch_threads"])
    torch.manual_seed(1337)

    config.ensure_dirs()
    print("NWP Emulator training")
    print(config.describe())
    print(f"  torch          : {torch.__version__}")
    print(resources.describe(RESOURCE_PLAN))

    # --- Data --------------------------------------------------------------
    run_folders = sorted([d for d in config.TENSOR_DIR.iterdir() if d.is_dir()]) \
        if config.TENSOR_DIR.exists() else []
    if not run_folders:
        run_folders = [config.TENSOR_DIR]      # flat layout fallback

    dataset = AutoregressiveDataset3D(run_folders)
    if len(dataset) == 0:
        print("\nNo training pairs found. Run the ingestion step first.")
        return 1

    n_val = max(1, int(len(dataset) * args.val_split)) if len(dataset) > 4 else 0
    n_train = len(dataset) - n_val
    if n_val:
        train_set, val_set = random_split(
            dataset, [n_train, n_val],
            generator=torch.Generator().manual_seed(1337))
    else:
        train_set, val_set = dataset, None

    train_loader = DataLoader(train_set, batch_size=args.batch_size, shuffle=True,
                              num_workers=RESOURCE_PLAN["dataloader_workers"])
    val_loader = DataLoader(val_set, batch_size=args.batch_size,
                            num_workers=RESOURCE_PLAN["dataloader_workers"]) \
        if val_set else None

    print(f"  train / val    : {n_train} / {n_val} transitions")

    # --- Model -------------------------------------------------------------
    model = NWPEmulator3D()
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    weights = torch.tensor(config.LOSS_WEIGHTS)

    n_params = sum(p.numel() for p in model.parameters())
    print(f"  parameters     : {n_params:,}")

    start_epoch, best_loss = 0, float("inf")
    latest = config.CHECKPOINT_DIR / "latest.pth"
    if args.resume and latest.exists():
        ckpt = torch.load(latest, weights_only=False)
        model.load_state_dict(ckpt["model_state"])
        optimizer.load_state_dict(ckpt["optimizer_state"])
        start_epoch = ckpt["epoch"] + 1
        best_loss = ckpt.get("best_loss", float("inf"))
        print(f"  resumed        : from epoch {start_epoch}")

    if args.dry_run:
        state_t, state_t1 = next(iter(train_loader))
        print(f"\nDry run -- input {tuple(state_t.shape)} -> target {tuple(state_t1.shape)}")
        out = model(state_t)
        loss, ch = control_variable_loss(out, state_t1, weights)
        print(f"  output {tuple(out.shape)} | loss {loss.item():.6f}")
        print("  per-channel:", {c: round(v.item(), 5) for c, v in zip(config.CHANNELS, ch)})
        print("\nWiring verified.")
        return 0

    # --- Train -------------------------------------------------------------
    print()
    for epoch in range(start_epoch, args.epochs):
        t0 = time.time()
        train_loss, train_ch = run_epoch(model, train_loader, optimizer, weights, train=True)

        msg = f"epoch {epoch + 1:3d}/{args.epochs}  train {train_loss:.6f}"
        if val_loader:
            val_loss, _ = run_epoch(model, val_loader, optimizer, weights, train=False)
            msg += f"  val {val_loss:.6f}"
            track = val_loss
        else:
            track = train_loss
        msg += f"  ({time.time() - t0:.1f}s)"
        print(msg)
        print("   channels:", {c: round(v.item(), 5) for c, v in zip(config.CHANNELS, train_ch)})

        save_checkpoint(latest, model, optimizer, epoch, best_loss)
        if track < best_loss:
            best_loss = track
            save_checkpoint(config.CHECKPOINT_DIR / "best.pth",
                            model, optimizer, epoch, best_loss)
            print(f"   new best ({best_loss:.6f}) -> best.pth")

    print(f"\nDone. Best loss {best_loss:.6f}. Checkpoints in {config.CHECKPOINT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
