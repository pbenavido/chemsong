import pandas as pd
import pretty_midi

# Example DataFrame with bond energies
data = {
    "Bond Type": ["C-C (aromatic)", "C-C (aliphatic)", "C-H", "C-Cl"],
    "Energy (kJ/mol)": [615, 348, 414, 328],
}
df = pd.DataFrame(data)

# Function to map energy values to MIDI notes
def energy_to_midi_note(energy, scale):
    # Simple linear mapping, can be more complex based on your data
    return scale[int(energy % len(scale))]


# Define a C Arabic minor scale (Hijaz)
c_arabic_minor_scale = [
    60,
    61,
    64,
    65,
    67,
    68,
    71,
]  # MIDI notes for C4, D♭4, E4, F4, G4, A♭4, B4

# Create a PrettyMIDI object
music = pretty_midi.PrettyMIDI()

# Create an instrument instance for a Piano (or any other instrument)
piano = pretty_midi.Instrument(
    program=pretty_midi.instrument_name_to_program("Acoustic Grand Piano")
)

# Add notes to the instrument
time = 0  # Start at the beginning
for index, row in df.iterrows():
    note_number = energy_to_midi_note(row["Energy (kJ/mol)"], c_arabic_minor_scale)
    note = pretty_midi.Note(velocity=100, pitch=note_number, start=time, end=time + 1)
    piano.notes.append(note)
    time += 1

# Add the instrument to the PrettyMIDI object
music.instruments.append(piano)

# Save the MIDI file
music.write("chemical_bonds_hijaz.mid")
