#!/usr/bin/env python3

import sys

COMMANDS = {
    "binding-frequency": "martinisoup.scripts.binding_frequency",
    "residence-times": "martinisoup.scripts.residence_times",
    "msd": "martinisoup.scripts.msd",
    "msd-fitter": "martinisoup.scripts.msd_fitter",
}

def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print("Usage: martinisoup <command> [args...]")
        print("\nCommands:")
        for name in COMMANDS:
            print(f"  {name}")
        sys.exit(0 if len(sys.argv) >= 2 else 1)

    command = sys.argv[1]
    if command not in COMMANDS:
        print(f"Unknown command: {command}")
        print(f"Available commands: {', '.join(COMMANDS)}")
        sys.exit(1)

    # Rewrite argv so the subcommand's argparse sees the right program name
    sys.argv = [f"martinisoup {command}"] + sys.argv[2:]

    import importlib
    mod = importlib.import_module(COMMANDS[command])
    mod.main()
