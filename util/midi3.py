import time

import mido
import pandas as pd
import pyo

from util.bond_energies import bond_energies
from util.scale_reference import *

# Initialize Pyo server
s = pyo.Server().boot()
s.start()

# Data Processing Function
def normalize_data(bond_df: pd.DataFrame) -> pd.DataFrame:

    sorted_bond_energies = sorted(set(bond_energies.values()))

    # normalize energies for midi mapping
    bond_df["Energy Index"] = bond_df[
        "Bond Energies"
    ].apply(  # TODO: change to "make sure key ["energies"] matches the key in the dataframe"
        lambda energies: [sorted_bond_energies.index(e) for e in energies]
    )

    return bond_df


# MIDI Generation Function
def generate_midi(data: list[int]) -> list[mido.Message]:
    midi_messages = []
    for value in data:
        note = map_to_scale(value)
        midi_messages.append(mido.Message("note_on", note=note, velocity=64, time=0))
        midi_messages.append(mido.Message("note_off", note=note, velocity=64, time=500))
    return midi_messages


# Map Numbers to C Minor Scale
def map_to_scale(note: int) -> int:
    # Mapping the number to a note in scale selected by tkinter dropdown

    scale = a_minor()

    # get dropdown value from where it's written:
    with open("util/dropdown_value.txt", "r") as f:
        scale_selection = f.read()
        f.close()

    if scale_selection == "a_major":
        scale = a_major()
    elif scale_selection == "a_minor":
        scale = a_minor()
    elif scale_selection == "a_minor_h":
        scale = a_minor_h()
    elif scale_selection == "b_major":
        scale = b_major()
    elif scale_selection == "b_minor":
        scale = b_minor()

    return scale[note % len(scale)]


# Synthesize Sine Wave with Reverb
# play with this function to get the sound you want
def synthesize_sound(midi_note: mido.Message):
    freq = midi_to_frequency(midi_note.note)
    osc = pyo.Sine(freq, mul=0.5).out()
    reverb = pyo.Freeverb(osc, size=0.5, damp=0.5, bal=0.3).out()

    # Play the note for a specific duration, then stop
    play_duration = 0.25  # Duration in seconds
    time.sleep(play_duration)

    # Stop the sound
    osc.stop()
    reverb.stop()


# Main synth function
def df_to_notes(df: pd.DataFrame):
    # normalize data
    normalized_dataframe = normalize_data(df)

    # play all notes in dataframe, iterating through steps
    for step in normalized_dataframe["Energy Index"]:
        time.sleep(0.25)  # delay between steps
        midi_messages = generate_midi(step)
        for note in midi_messages:
            # Send MIDI messages to the synthesizer
            synthesize_sound(note)
