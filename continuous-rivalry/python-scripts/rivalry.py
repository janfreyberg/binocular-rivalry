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

stimsize = 8  # size of stimuli in deg of visual angle
repetitions = 2  # how many trials per condition

screenrefresh = 60.0  # refresh rate of screen
frequencies = [5.0, 7.5]  # the two (competing) frequencies

simulationpercepts = [1.5, 3.5]  # interval for simulated percept length
# durations will be sampled evenly from this interval

wedges = [8, 16]  # how many "pie" slices in stimulus
concentrics = [6, 3]  # how many "rings" in stimulus

colors = [1, 0, 0.9]  # the rgb layers of your stimulus
# the 0.8 is due to the anaglyph blue being brighter than red

# What buttons you want to use, and what image they correspond to:
responses = {'right': 'red', 'up': 'mixture', 'left': 'blue'}

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
square = visual.Circle(win, radius=np.sqrt(2) * stimsize / 2,
                       fillColor=None, autoDraw=True, edges=128,
                       lineWidth=2, pos=[0, 0], lineColor=-1)
# Make fixation cross (using list comprehension)
fixation = visual.GratingStim(win, tex="sqr", mask="cross", sf=0, size=0.5,
                              pos=[0, 0], color=-1,
                              autoDraw=True)
# Make a dummy message
message = visual.TextStim(win, units='norm', pos=[0, 0.75], height=0.06,
                          alignVert='center', alignHoriz='center',
                          text='')


# The arrays needed to construct the stimuli:
# Mask
mask = np.concatenate([np.array([-1]), np.ones(min(concentrics))])
# Checkerboards
texlayer = [np.repeat(
    np.repeat(
            np.tile(np.array([[1, 0], [0, 1]]), (concentrics[stim] +
                                                 concentrics[stim] /
                                                 min(concentrics),
                                                 wedges[stim] / 2)),
            concentrics[stim - 1] + concentrics[stim - 1] /
            min(concentrics), axis=0
            ),
    6 * wedges[stim - 1], axis=1
) * color * 2 - 1
    for stim, color in enumerate([colors[0], colors[2]])]
# Blank (all black)
blanklayer = np.zeros(texlayer[0].shape) - 1

# Make stimuli. This will be remade during each trial anyways tho
stimuli = visual.RadialStim(win, size=stimsize, pos=fixations[stim].pos,
                            mask=mask, colorSpace='rgb', angularRes=600,
                            # Radial and Angular Frequencies:
                            radialCycles=0.5,
                            angularCycles=1,
                            # Baseline Texture:
                            tex=np.dstack((texlayer[0],
                                           blanklayer,
                                           texlayer[1]))
                            )


# Define an instruction function
def instruct(displaystring):
    message.text = displaystring + "\nPress [space] to continue."
    message.draw()
    win.flip()
    keyboard.reporting = True
    if 'escape' in keyboard.waitForReleases(keys=[' ', 'escape']):
        raise KeyboardInterrupt
    keyboard.reporting = False


def trigger(value):
    print(value)


# function for demonstration
def demonstrate():

    core.wait(0.5)
    # Draw Red:
    message.text = ("This is what a red stimulus will look like. In "
                    "the experiment, you would report seeing this by pressing "
                    "right.\nPress [space] to continue.")
    frame = 0
    keyboard.reporting = True
    while ' ' not in keyboard.getReleases(keys=[' ']):
        frame += 1  # update frame tally
        contrmod = [np.sign(np.cos(np.pi * frame / (screenrefresh /
                                                    frequencies[layer])))
                    for layer in range(2)]
        stimuli.tex = np.dstack((texlayer[0] * contrmod[0],
                                 blanklayer,
                                 blanklayer))
        message.draw()
        stimuli.draw()
        win.flip()
    # Clear & Flush
    win.flip()
    keyboard.getKeys()
    keyboard.reporting = False
    core.wait(1)

    # Draw Blue:
    message.text = ("And this is what a blue stimulus will look like. "
                    "In the experiment, you would report seeing this by "
                    "pressing left.\nPress [space] to continue.")
    frame = 0
    keyboard.reporting = True
    while ' ' not in keyboard.getReleases(keys=[' ']):
        frame += 1  # update frame tally
        contrmod = [np.sign(np.cos(np.pi * frame / (screenrefresh /
                                                    frequencies[layer])))
                    for layer in range(2)]
        stimuli.tex = np.dstack((blanklayer,
                                 blanklayer,
                                 texlayer[1] * contrmod[1]))
        message.draw()
        stimuli.draw()
        win.flip()
    # Clear & Flush
    win.flip()
    keyboard.getKeys()
    keyboard.reporting = False
    core.wait(1)

    # Draw both rotating:
    message.text = ("And this is what the two will look like together. "
                    "In the experiment, you will have to report what "
                    "you are seeing at any given point by pressing "
                    "either [left], [right], or [up]."
                    "\nPress [space] to continue.")
    frame = 0
    keyboard.reporting = True
    while not keyboard.getReleases(keys=[' ']):
        frame += 1  # update frame tally
        contrmod = [np.sign(np.cos(np.pi * frame / (screenrefresh /
                                                    frequencies[layer])))
                    for layer in range(2)]
        stimuli.tex = np.dstack((texlayer[0] * contrmod[0],
                                 blanklayer,
                                 texlayer[1] * contrmod[1]))
        message.draw()
        stimuli.draw()
        win.flip()
    # Clear & Flush
    win.flip()
    keyboard.reporting = False
    core.wait(1)


# Define a trial function:
def rivaltrial(info):

    # Arrays to keep track of responses:
    keyArray = np.zeros((int(trialdur * screenrefresh), 3),
                        dtype=int)  # start with up pressed
    timeArray = np.zeros((int(trialdur * screenrefresh), 1),
                         dtype=float)  # start at time 0

    # Wait for participant to start:
    message.text = ("Trial number " + str(1 + trials.thisN) + " of " +
                    str(trials.nTotal) +
                    ". Remember: UP is for mixture, LEFT is for "
                    "blue, and RIGHT is for red.\nPress and hold [up] "
                    "to begin the trial.")
    message.draw()
    win.flip()
    keyboard.reporting = True
    if "escape" in keyboard.waitForPresses(keys=['up', 'escape']):
        raise KeyboardInterrupt("You interrupted the script manually.")

    # This will keep track of response time:
    trialClock = core.Clock()
    win.callOnFlip(trialClock.reset)

    # Make a list of contrast modulators to save time during trial
    contrmod = [[(int(np.sign(np.cos(
        np.pi * frame / (screenrefresh / info['frequency'][layer])))) + 1) / 2
        for layer in range(2)]
        for frame in range(int(trialdur * screenrefresh))]

    # If this is a simulation trial, make a list to determine stimulus
    stimsequence = np.zeros((int(trialdur * screenrefresh), 1),
                            dtype=int) + 2
    if info['simulation']:
        perceptstart = 0
        currentpercept = 0
        while perceptstart < stimsequence.size:
            perceptframes = int(screenrefresh * (np.random.random() *
                                                 (simulationpercepts[1] -
                                                  simulationpercepts[0]) +
                                                 simulationpercepts[1]))
            # assign the percept
            if perceptstart + perceptframes < stimsequence.size:
                stimsequence[
                    perceptstart:(perceptstart + perceptframes)
                ] = currentpercept
                # update next percept start and value
                currentpercept = (currentpercept + 1) % 2
                perceptstart += perceptframes
            else:
                stimsequence[
                    perceptstart:(stimsequence.size)
                ] = currentpercept
                perceptstart = stimsequence.size

    # Pre-calculate all possible textures
    textures = [
        [[np.dstack((blanklayer,
                     blanklayer,
                     texlayer[info['pattern'][1]] *
                     bluecontrast))
          for redcontrast in [-1, 1]]
         for bluecontrast in [-1, 1]],
        [[np.dstack((texlayer[info['pattern'][0]] *
                     redcontrast,
                     blanklayer,
                     blanklayer))
          for redcontrast in [-1, 1]]
            for bluecontrast in [-1, 1]],
        [[np.dstack((texlayer[info['pattern'][0]] *
                     redcontrast,
                     blanklayer,
                     texlayer[info['pattern'][1]] *
                     bluecontrast))
          for redcontrast in [-1, 1]]
         for bluecontrast in [-1, 1]]
    ]

    # Make the textures
    stimarray = [[[visual.RadialStim(win, size=stimsize,
                                     pos=fixations[0].pos,
                                     mask=mask, colorSpace='rgb',
                                     angularRes=600,
                                     # Radial and Angular Frequencies:
                                     radialCycles=0.5,
                                     angularCycles=1,
                                     # Baseline Texture:
                                     tex=textures[stim][red][blue]
                                     )
                   for red in range(2)]
                  for blue in range(2)]
                 for stim in range(3)]

    # Now cycle through frames (actual trial):
    trigger(100 + trials.thisN)
    for frame in range(int(trialdur * screenrefresh)):

        # Update the stimulus texture
        stimarray[
            stimsequence[frame]
        ][contrmod[frame][0]][contrmod[frame][1]].draw()
        # Flip
        win.flip()

        # Check response & save to array
        keyArray[frame, :] = [response in keyboard.state.keys()
                              for response in responses.keys()]
        timeArray[frame] = trialClock.getTime()

        # if response is different, trigger EEG
        if not np.array_equal(keyArray[frame, :], keyArray[frame - 1, :]):
            # send value as decimal rep of binary vector
            trigger(1 + sum(1 << i
                            for i, b in enumerate(keyArray[frame, :])
                            if b))

    # After trial is over, clear the screen
    win.flip()
    keyboard.reporting = False

    # Save the trial info:
    w = csv.DictWriter(open(datadir +
                            'trial-%03d-info.csv' %
                            trials.thisN, 'w+'),
                       info.keys())
    w.writeheader()
    w.writerow(info)
    # Save the data:
    np.savetxt(datadir + 'trial-%03d-data.csv' % trials.thisN,
               np.hstack([timeArray, keyArray, stimsequence]), fmt='%3.5f',
               delimiter=', ',
               header=', '.join(['time'] +
                                responses.values() +
                                ['simulationsequence']))

    # return the two arrays
    return keyArray, timeArray


# Define a break function:
def rivalbreak():
    fixation.autoDraw = False  # disable fixation
    # Display message
    keyboard.reporting = True
    for sec in range(int(breakdur)):
        message.text = "Break for %d seconds." % (breakdur - sec)
        message.draw()
        win.flip()
        # Pause for a second (interrupt if key pressed)
        if 'escape' in keyboard.waitForReleases(maxWait=1, keys=['escape']):
            # Escape stops script:
            raise KeyboardInterrupt("You interrupted the script manually.")
    keyboard.reporting = False
    fixation.autoDraw = True  # enable fixation
    win.flip()


# Explain the experiment
instruct("This experiment is called binocular rivalry, and you will need the "
         "red and blue glasses to do so. Your right eye should be viewing "
         "through a red lense, and your left eye should be viewing through "
         "a blue lens.")
instruct("In this experiment, we will be showing a different pattern "
         "to each of your eyes. This will mean that for a while, you will "
         "be seeing one image, and then it will suddenly switch.")
instruct("All you need to do during each trial of this experiment is "
         "report what you are seeing continuously. To do so, you have "
         "to hold down either the UP arrow, the LEFT arrow, or the RIGHT "
         "arrow.")
instruct("These buttons stand for different images:\nThe RIGHT arrow stands "
         "for the RED image.\nThe LEFT arrow stands for the blue image.\n"
         "And the UP arrow stands for a mixture, so for example, when about "
         "50% is red and 50% is blue.")
instruct("It's rare for an image to be 100% blue or red, but you should "
         "still only press the UP arrow when you can't decide which color "
         "is dominant.")
instruct("Please report what you are seeing continuously. So, during one trial,"
         " you should be pressing one button at any one time. Next, we will "
         "demonstrate what the images are like.")

# demonstrate the stimuli:
demonstrate()

# check before starting:
instruct("If you have any other questions, please ask the experimenter now. "
         "Otherwise, go ahead to start the experiment.")

# randomise trial sequence:
trials = data.TrialHandler(
    [{'pattern': [0, 1] if pattern else [1, 0],
      'frequency': frequencies if flicker else list(reversed(frequencies)),
      'flickering': False if trialtype == 0 else True,
      'simulation': True if trialtype == 2 else False}
     for pattern in range(2)
     for flicker in range(2)
     for trialtype in range(3)],
    repetitions)

# make a directory for data storage
datadir = os.path.join(os.getcwd(), 'data',  sessionInfo['time'] +
                       ' ' + sessionInfo['subject'], '')
os.makedirs(datadir)

# Loop through trials:
for trial in trials:
    rivaltrial(trial)
    # Take a break (not on last)
    if trials.thisN < trials.nTotal:
        rivalbreak()

# Close everything neatly
win.close()
tracker.destroy()
core.quit()
