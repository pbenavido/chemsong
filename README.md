# Chemsong

## Overview

Chemsong is a program that combines chemical reaction visualization with audio mappings. The code translates the energies of chemical bonds into musical scales so that chemical processes can be observed along with a related audio sequence.

## Features

- **Chemical Reaction Visualization**: Uses RDKit to render molecular structures and transformations during chemical reactions.
- **Bond Energy to Audio Mapping**: Translates chemical bond energies into notes on musical scales, creating an auditory representation.
- **Data Analysis and Visualization**: Employs Pandas for data handling and RDKit.Draw for mapping chemical properties to musical notes.

## Tools

- **RDKit**: Chosen for its molecular visualization capabilities and ease of mapping chemical properties.

## Project Structure

- `app/`: central scripts for running Chemsong
- `util`: 'under-the-hood' Chemsong subsystems .
- `example/`: Scripts showcasing RDKit visualizations and PrettyMIDI audio mappings. Generally speaking, these scripts don't work.
- `experiment/`: Initial explorations with RDKit, ChemTools, and ChemPy. Some of these work I think.

## Future Development

- **Refine Audio Mapping**: Enhance the translation of bond energies to musical scales.
- **Expand Reaction Range**: Include a wider variety of chemical reactions.
- **Interactive Features**: Develop interactive elements for real-time manipulation and audio feedback.

## Getting Started

Requirements include Python (3.7+ recommended), RDKit, Pandas, and related libraries used for audio processing. See `requirements.txt` for details.

## Contributing

Contributions are welcome! Enhance audio mappings, visualizations, or add new features to enrich the project.

## License

Open-sourced under the [MIT License](LICENSE).
