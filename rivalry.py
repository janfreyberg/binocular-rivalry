from psychopy import core, visual, event, monitors, gui, data
from psychopy.iohub import launchHubServer  # for parallel key collection
from datetime import datetime
import random  # random for randomisation
import numpy as np  # for array handling
import os
import os.path
import csv

# Stimulus Parameters
trialdur = 5.0  # trial duration in seconds
breakdur = 5.0  # minimum duration of a break in seconds
contrast = 0.6  # the contrast of the gratings
stimsize = 4  # size of stimuli in deg of visual angle
spatialfrequency = 2.0
# how many trials per condition
# nb multiply this by 8 to get trialnum
repetitions = 1

# Checkerboard stimuli
screenrefresh = 60.0
frequencies = [5.0, 7.5]
wedges = [8, 16]
concentrics = [6, 3]
colors = [(0.8, 0, 0), (0, 0.8, 0)]
# calculate how much is masked in the middle
maskedcircles = np.divide(concentrics, min(concentrics))

# What buttons you want to use, and what image they correspond to:
responses = {'right': 'red', 'up': 'mixture', 'left': 'green'}

# Get participant information:
sessionInfo = {'subject': 'test',
               'time': datetime.now().strftime('%Y-%m-%d %H-%M')}
assert gui.DlgFromDict(dictionary=sessionInfo, title='Pupil Dilation').OK

# Set the screen parameters: (This is important!)
screen = monitors.Monitor("testMonitor")
screen.setSizePix([1680, 1050])
screen.setWidth(47.475)
screen.setDistance(57)
# Open the display window:
win = visual.Window([500, 500], allowGUI=True, monitor=screen,
                    units='deg', fullscr=True)

# Start the keyboard monitor device
io = launchHubServer()
keyboard = io.devices.keyboard
keyboard.reporting = True

# Make the vergence cues (using list comprehension)
squares = [visual.Circle(win, radius=np.sqrt(2) * stimsize / 2,
                         fillColor=None, autoDraw=True, edges=128,
                         lineWidth=2, pos=[side * stimsize, 0], lineColor=-1)
           for side in [-1, 1]]
# Make fixation cross (using list comprehension)
fixations = [visual.GratingStim(win, tex="sqr", mask="cross", sf=0, size=0.5,
                                pos=[side * stimsize, 0], color=-1,
                                autoDraw=True)
             for side in [-1, 1]]


# The arrays needed to construct the stimuli:
mask = np.concatenate([np.array([-1]), np.ones(min(concentrics))])
relevantcolor = np.array([[1, 1, -1, -1], [-1, -1, 1, 1],
                          [-1, -1, 1, 1], [1, 1, -1, -1]])


def trigger(value):
    print(value)


def findVergence():
    # First make the vergence cues distinct:
    fixations[0].color = squares[0].lineColor = 1
    win.flip()
    # These are the distances we'll move the stimuli by:
    distances = {'up': [0.2, 0], 'down': [-0.2, 0],
                 'right': [0.05, 0], 'left': [-0.05, 0],
                 'space': [stimsize / 2, 0], ' ': [stimsize / 2, 0]}
    # Now use button presses to move them:
    while True:
        keys = keyboard.waitForPresses(keys=distances.keys())
        if 'space' in keys or ' ' in keys:
            # make color the same
            fixations[0].color = squares[0].lineColor = -1
        elif 'return' in keys or 'end' in keys:
            # exit loop
            return
        elif 'escape' in keys:
            raise KeyboardInterrupt("You interrupted the script manually.")
        # Update the stimulus positions:
        squares[0].pos = fixations[0].pos = \
            fixations[0].pos - distances[keys[0].key]
        squares[1].pos = fixations[1].pos = \
            fixations[1].pos + distances[keys[0].key]
        win.flip()


# function for demonstration
def demonstrate():
    # Make some gratings:
    gratings = [visual.RadialStim(win, size=stimsize, pos=fixations[stim].pos,
                                  mask=mask, colorSpace='rgb', angularRes=600,
                                  # Radial and Angular Frequencies:
                                  radialCycles=(concentrics[stim] +
                                                concentrics[stim] /
                                                min(concentrics)) / 2,
                                  angularCycles=wedges[stim] / 2,
                                  # Baseline Texture:
                                  tex=np.stack([(relevantcolor + 1) * rgb - 1
                                                for rgb in colors[stim]],
                                               axis=2))
                for stim in range(2)]

    # Draw the same grating both sides:
    for stim in range(2):
        gratings[stim].draw()  # draw on side one
        # swap location:
        gratings[0].pos, gratings[1].pos = gratings[1].pos, gratings[0].pos
        gratings[stim].draw()  # draw on side two
        win.flip()
        keyboard.waitForPresses()

    # Draw the gratings rivalling:
    [grating.draw() for grating in gratings]
    win.flip()
    keyboard.waitForPresses()
    # clear the screen:
    win.flip()


# Define a trial function:
def rivaltrial(info):
    # Make gratings for this trial:
    gratings = [visual.RadialStim(win, size=stimsize, pos=fixations[stim].pos,
                                  mask=mask, colorSpace='rgb', angularRes=600,
                                  # Radial and Angular Frequencies:
                                  radialCycles=(info['concentrics'][stim] +
                                                info['concentrics'][stim] /
                                                min(concentrics)) / 2,
                                  angularCycles=info['wedges'][stim] / 2,
                                  # Baseline Texture:
                                  tex=np.stack([(relevantcolor + 1) * rgb - 1
                                                for rgb in info['color'][stim]],
                                               axis=2))
                for stim in range(2)]

    # Arrays to keep track of responses:
    keyArray = np.array([[0.0, 1.0, 0.0]])
    timeArray = np.array([[0.0]])

    # Wait for participant to start:
    if "escape" in keyboard.waitForPresses(keys=['up', 'escape']):
        raise KeyboardInterrupt("You interrupted the script manually.")

    # This will keep track of response time:
    trialClock = core.Clock()
    win.callOnFlip(trialClock.reset)

    # Now cycle through frames:
    for frame in range(int(trialdur * screenrefresh)):

        # Update phase
        [gratings[stim].setAngularPhase(
            0.5 * (info['flickering'] and
                   frame % (screenrefresh / (2 * info['frequency'][stim])) == 0)
        ) for stim in range(2)]

        # Draw and flip
        [grating.draw() for grating in gratings]
        win.flip()

        # Check response
        keyArray = np.vstack((keyArray,
                              np.array([[response in keyboard.state.keys()
                                         for response in responses.keys()]])))
        timeArray = np.vstack([timeArray, [[trialClock.getTime()]]])

        # if response is different, trigger
        if not np.array_equal(keyArray[-1, :], keyArray[-2, :]):
            # send value as decimal rep of binary vector
            trigger(1 + sum(1 << i for i, b in enumerate(keyArray[-1, :]) if b))

        if "escape" in keyboard.state.keys():
            raise KeyboardInterrupt("You interrupted the script manually.")

    win.flip()  # clear the screen

    # Save the trial info:
    w = csv.DictWriter(open(datadir +
                            'trial-%03d-info.csv' %
                            trials.thisN, 'w+'),
                       info.keys())
    w.writeheader()
    w.writerow(info)
    # Save the data:
    np.savetxt(datadir + 'trial-%03d-keys.csv' % trials.thisN,
               np.hstack([timeArray, keyArray]), fmt='%3.5f',
               delimiter=', ',
               header=', '.join(['time'] +
                                responses.values()))

    # return the two arrays
    return keyArray, timeArray


# Define a break function:
def rivalbreak():
    fixations[0].autoDraw = fixations[1].autoDraw = False  # disable fixation
    # Display message
    breakmsgs = [visual.TextStim(win, units='deg', wrapWidth=stimsize,
                                 alignVert='center', alignHoriz='center',
                                 height=0.5, color=-1, text="",
                                 pos=fixations[side].pos - [0, 0.25 * stimsize])
                 for side in range(2)]
    for sec in range(int(breakdur)):
        breakmsgs[0].text = breakmsgs[1].text = ("Break for %d seconds."
                                                 % (breakdur - sec))
        [breakmsg.draw() for breakmsg in breakmsgs]
        win.flip()
        # Pause for a second (interrupt if key pressed)
        breakkeys = keyboard.waitForPresses(maxWait=1,
                                            keys=['escape', 'return'])
        if breakkeys and 'escape' in breakkeys:
            # Escape stops script:
            raise KeyboardInterrupt("You interrupted the script manually.")
        elif breakkeys and 'return' in breakkeys:
            break  # Enter skips the break:
    breakmsgs[0].text = breakmsgs[1].text = ("Begin the next trial "
                                             "with up key.")
    [breakmsg.draw() for breakmsg in breakmsgs]
    win.flip()
    fixations[0].autoDraw = fixations[1].autoDraw = True  # enable fixation


# find the point at which the stimuli merge:
# findVergence()
# demonstrate the stimuli:
# demonstrate()

# randomise trial sequence:
trials = data.TrialHandler(
    [{'color': colors if color else list(reversed(colors)),
      'wedges': wedges if pattern else list(reversed(wedges)),
      'concentrics': concentrics if pattern else list(reversed(concentrics)),
      'frequency': frequencies if flicker else list(reversed(frequencies)),
      'flickering': False if flicker == 2 else True,
      'simulation': True if flicker == 2 else True}
     for pattern in range(2)
     for color in range(2)
     for flicker in range(4)],
    repetitions)

print(trials.nTotal)

# make a directory for data storage
datadir = os.path.join(os.getcwd(), 'data',  sessionInfo['time'] +
                       ' ' + sessionInfo['subject'], '')
os.makedirs(datadir)

# Loop through trials:
keyArrays = []
timeArrays = []
for trial in trials:
    (keyArray, timeArray) = rivaltrial(trial)
    # Take a break (not on last)
    if trials.thisN < trials.nTotal:
        rivalbreak()
