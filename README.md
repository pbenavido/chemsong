# Chemical Reaction Visualization with Audio Mapping

## Project Overview

This project explores visualizing chemical reactions and maps the energies of chemical bonds to musical scales, creating an auditory experience of chemical processes. 

### Key Features

- **Visualization of Chemical Reactions**: Utilize advanced chemical informatics tools to visualize molecular structures and their transformations during chemical reactions.
- **Audio Mapping of Bond Energies**: Convert the energies of chemical bonds into corresponding notes on a musical scale, providing an auditory representation of chemical bonds.
- **Data Analysis with Pandas and RDKit**: Leverage Pandas dataframes and RDKit.Draw for efficient data handling and mapping between chemical properties and musical notes.

## Tools Explored

During the initial phase of the project, various tools were explored to determine their suitability:

- **RDKit**: Emerged as the primary tool for this project due to its robust capabilities in molecular visualization and ease of mapping chemical properties to data structures.
- **ChemTools and ChemPy**: Initially considered for their specialized capabilities in quantum chemistry and chemical kinetics, but found to be more oriented towards stoichiometry and electronic structure analysis, which are less directly applicable to the project's goals.

## Choosing RDKit and PrettyMIDI

After exploring various options, RDKit and PrettyMIDI were chosen as the primary tools due to:

- **Direct Bond Analysis**: RDKit's ability to analyze and categorize bond types directly.
- **Integration with Pandas**: Good integration with Pandas dataframes, facilitating easy manipulation and mapping of chemical data.
- **Visualization Strengths**: Great support for molecular structure visualization, crucial for this project's core objective.
- **Audio Mapping with PrettyMIDI**: PrettyMIDI's flexibility in creating and manipulating MIDI files, allowing for the mapping of chemical properties to musical notes.

## Project Structure

- **`example/` Directory**: Contains functional scripts demonstrating the core features of the project, including RDKit visualizations and PrettyMIDI audio mappings. These serve as 'unit demos' and will be expanded into a full application in the `app/` directory.
- **`experiment/` Directory**: Showcases initial explorations and demonstrations of various chemistry tools, including RDKit, ChemTools, and ChemPy.

## Future Directions

- **Enhancing Audio Mapping**: Further refine the mapping of bond energies to musical scales using PrettyMIDI to create more nuanced and varied auditory representations.
- **Expanding Chemical Reactions**: Broaden the range of chemical reactions visualized and sonified, exploring both organic and inorganic chemistry.
- **Interactive Features**: Develop interactive features that allow users to manipulate molecular structures and hear corresponding changes in the music.

## Getting Started

To get started with this project, you'll need to have Python installed (3.7+ recommended), along with RDKit, Pandas, PrettyMIDI, and any necessary audio processing libraries (see `requirements.txt`). Detailed setup instructions will be provided in subsequent sections.

## Contribution

Contributions to this project are welcome! Whether it's adding new features, refining the audio mapping algorithms, or enhancing the visualizations, your input can help take this project to new heights.

## License

This project is open-sourced under the [MIT License](LICENSE).
