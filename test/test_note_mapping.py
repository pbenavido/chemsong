from ..util.midi3 import normalize_data, df_to_notes

from loguru import logger
import pandas as pd
test_df = pd.DataFrame({
    'Step': [1, 2],
    'Bond Energies': [[835], [346, 835]]
})

def test_note_mapping():
    normalized_data = normalize_data(test_df)
    logger.info(f"\n{normalized_data}")
    
    assert normalized_data["Energy Index"][:2][0] == [47]
    assert normalized_data["Energy Index"][:2][1] == [28, 47]

def test_note_playing():
    df_to_notes(test_df)
    assert True
