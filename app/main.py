from util.reaction import molecules_to_bond_energy_df
from util.visualize import render_smiles, verify_smiles
from util.midi import df_to_notes
import tkinter as tk
from tkinter import scrolledtext
from tkinter import messagebox
from PIL import Image, ImageTk
from rdkit import Chem

def chemsong(mols):
    # get bond energies from mols
    bond_df = molecules_to_bond_energy_df(mols)

    # visualize steps and display in Tkinter window
    image_labels = []  # List to store image labels
    for step, smiles in mols.items():
        img = render_smiles(smiles)
        photo = ImageTk.PhotoImage(img)
        label = tk.Label(image_frame, image=photo)
        label.image = photo  # Keep a reference
        label.pack()
        image_labels.append(label)  # Store the label

    # play notes
    df_to_notes(bond_df)


def process_reaction(entries):
    mols = {}
    for i, entry in enumerate(entries):
        step = i + 1
        smiles = entry.get()
        if not verify_smiles(smiles):
            messagebox.showerror("Error", f"Invalid SMILES notation at step {step}")
            return
        # Split the smiles string into a list and assign to the step
        mols[step] = smiles.split()
    chemsong(mols)


def add_step():
    step_number = len(steps_entries) + 1
    entry = tk.Entry(entry_frame, width=30)  # Adjust the width of the entry
    entry.grid(row=step_number+1, column=0)
    steps_entries.append(entry)

# Read the SMILES guide from the text file
with open("smiles_notation_guide.txt", "r") as file:
    smiles_guide = file.read()

root = tk.Tk()
root.title("Chemsong")
root.geometry("800x600")  # Set a larger size for the window

# Title and Image at the top
top_frame = tk.Frame(root)
top_frame.pack(fill=tk.X)

title_label = tk.Label(top_frame, text="CHEMSONG", font=("Arial", 24))
title_label.pack(side=tk.LEFT, padx=10)

# Load and place the image
image = Image.open("chemsong.png")
image = image.resize((100, 100))  # Resize the image as needed
photo = ImageTk.PhotoImage(image)
image_label = tk.Label(top_frame, image=photo)
image_label.pack(side=tk.RIGHT, padx=10)

# Frame for entries and buttons
entry_frame = tk.Frame(root)
entry_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.Y)

steps_entries = []
tk.Button(entry_frame, text="Add Step", command=add_step).grid(row=0, column=0)
tk.Button(entry_frame, text="Process", command=lambda: process_reaction(steps_entries)).grid(row=1, column=0)

# Frame for SMILES guide
guide_frame = tk.Frame(root)
guide_frame.pack(side=tk.RIGHT, padx=10, pady=10, fill=tk.Y)
guide_text = scrolledtext.ScrolledText(guide_frame, wrap=tk.WORD, width=40, height=15)  # Adjust size
guide_text.insert(tk.INSERT, smiles_guide)
guide_text.config(state='disabled')
guide_text.pack()

# Frame for RDKit images
image_frame = tk.Frame(root)
image_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True, padx=10, pady=10)
root.mainloop()
