import pandas as pd

from util.bond_energies import bond_energies
from experiment.midi_draft_2 import map_to_a_minor_scale, play_notes_from_dataframe

# Example DataFrame
data = {
    0: [1,2,3],
    1: [[835],
        [346,835],
        [835,346]]
}
df = pd.DataFrame(data)

ordered_energy_list = list(bond_energies.values())

# normalize energy values using position in ordered_energy_list
for i in range(len(df)):
    for j in df[1][i]:
        norms = [ordered_energy_list.index(x) for x in df[1][i]]


# Play the notes
play_notes_from_dataframe(df)

# Note: This script won't produce sound in this environment. 
# You need to run it in your local environment with MIDI setup.
