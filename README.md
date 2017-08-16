# Binocular rivalry experiment scripts

This repository has some scripts for running binocular rivalry experiments.

You can run them in python (much nicer) or matlab (possibly has better timing...).

## Requirements

### Python requirements

You should be able to install everything required by running `pip install -r requirements.txt`.
If that doesn't work or you want to try using conda, just check out [its contents](requirements.txt) and install all mentioned there :)

### Matlab requirements

It's super hard to know which of matlabs various toolboxes functions come from. I would imagine the image processing toolbox is required.
You also need to [install psychtoolbox](http://psychtoolbox.org/).

You will also need the visual angle conversion scripts [here](https://github.com/autism-research-centre/convert-visangle).

## Different experiments

There are three experiments:

- `continuous-rivalry` is a standard ~40s rivalry experiment
- `trial-based-rivalry` is for testing binocular rivalry on lots of small (6s) trials
- `eeg-rivalry` is for flickering the stimuli to track them via the evoked frequency in EEG data
