import time

import mido
import pandas as pd


def map_to_a_minor_scale(note):
    # Mapping the number to a note in A minor scale
    a_minor_scale = [
        69,
        71,
        72,
        74,
        76,
        77,
        79,
    ]  # MIDI notes for A minor scale (A4 to G5)
    return a_minor_scale[note % len(a_minor_scale)]


def df_to_notes(df):
    # Setting up a MIDI port (this should be routed to a software synthesizer or DAW)
    with mido.open_output() as port:
        for index, row in df.iterrows():
            # Converting the numbers in the row to MIDI notes
            midi_notes = [map_to_a_minor_scale(note) for note in row]

            # Sending MIDI messages for each note
            for note in midi_notes:
                msg = mido.Message("note_on", note=note)
                port.send(msg)

            # Delay to let the notes play
            time.sleep(0.5)

            # Sending note_off messages to stop the notes
            for note in midi_notes:
                msg = mido.Message("note_off", note=note)
                port.send(msg)

            # Delay between rows
            time.sleep(0.5)
