"""
original values provided by Matt Fries, see moonpapasart.com/
"""

import math

# Frequencies list
FREAKwencies = [
    55.0,
    110.0,
    51.9131,
    103.826,
    61.7354,
    123.471,
    58.2705,
    116.541,
    32.7032,
    65.4064,
    36.7081,
    73.4162,
    34.6478,
    69.2957,
    41.2034,
    82.4069,
    38.8909,
    77.7817,
    43.6535,
    87.3071,
    48.9994,
    97.9989,
    46.2493,
    92.4986,
]

# Function to convert frequency to MIDI note number
def frequency_to_midi(frequency):
    freq = int(round(69 + 12 * math.log2(frequency / 440.0)))
    boosted_freq = freq + 30  # boosting to make audible
    return boosted_freq


# Function to convert MIDI note number to frequency
def midi_to_frequency(midi):
    return 440.0 * 2.0 ** ((midi - 69) / 12.0)


# Functions to return scales in MIDI note numbers
def a_minor():
    return [frequency_to_midi(FREAKwencies[i]) for i in [9, 11, 12, 14, 16, 17, 19, 21]]


def a_major():
    return [frequency_to_midi(FREAKwencies[i]) for i in [9, 11, 13, 14, 16, 18, 20, 21]]


def a_minor_h():
    return [frequency_to_midi(FREAKwencies[i]) for i in [9, 11, 12, 14, 16, 17, 20, 21]]


def b_minor():
    return [
        frequency_to_midi(FREAKwencies[i]) for i in [11, 13, 14, 16, 18, 19, 21, 23]
    ]


def b_major():
    return [
        frequency_to_midi(FREAKwencies[i]) for i in [11, 13, 15, 16, 18, 20, 22, 23]
    ]


scale_list = ["a_minor", "a_major", "a_minor_h", "b_minor", "b_major"]
