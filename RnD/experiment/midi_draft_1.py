import time

import numpy as np
import pygame
from loguru import logger

from util.bond_energies import bond_energies


def df_to_notes(bond_df):

    # Assuming bond_energies is already sorted
    sorted_bond_energies = sorted(bond_energies.items(), key=lambda item: item[1])
    sorted_bond_energies_dict = {k: i for i, (k, _) in enumerate(sorted_bond_energies)}

    # normalize energies for midi mapping
    bond_df["Normalized Energies"] = bond_df["Bond Energies"].apply(
        lambda energies: [sorted_bond_energies_dict.get(e) for e in energies]
    )

    # begin midi code section
    # Define the Arabic minor scale frequencies (in Hz)
    arabic_minor_scale = [
        261.63,
        293.66,
        311.13,
        349.23,
        392.00,
        415.30,
        466.16,
    ]  # C4 scale
    default_frequency = min(arabic_minor_scale)  # The lowest note in the scale

    def get_note_frequency(position, scale_length):
        if position is None:
            return default_frequency
        return arabic_minor_scale[position % scale_length]

    # Initialize pygame mixer
    pygame.mixer.init(frequency=44100, size=-16, channels=1)

    def generate_sound(frequency, duration=1):
        """Generate a sound of a specific frequency and duration."""
        sample_rate = 44100  # Samples per second
        t = np.linspace(
            0, duration, int(sample_rate * duration), endpoint=False
        )  # Time axis
        wave = np.sin(2 * np.pi * frequency * t)  # Sine wave
        buffer = (wave * 32767).astype(np.int16)  # Convert to 16-bit format
        return pygame.sndarray.make_sound(buffer)

    # Iterate over DataFrame and play notes
    for _, row in bond_df.iterrows():
        for pos in row["Normalized Energies"]:
            frequency = get_note_frequency(pos, len(arabic_minor_scale))
            sound = generate_sound(frequency)
            sound.play()
            time.sleep(0.5)  # Delay between notes

    # end midi code section
