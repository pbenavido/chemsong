import os
import random
import tkinter as tk
from tkinter import messagebox, scrolledtext, ttk

from PIL import Image, ImageTk
from rdkit import Chem
from ttkthemes import ThemedTk
from loguru import logger

from util.midi3 import df_to_notes
from util.reaction import molecules_to_bond_energy_df
from util.visualize import render_smiles, verify_smiles
from util.scale_reference import scale_list


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


def play_notes(entries):
    mols = {}
    for i, entry in enumerate(entries):
        step = i + 1
        smiles = entry.get()
        if not verify_smiles(smiles):
            messagebox.showerror("Error", f"Invalid SMILES notation at step {step}")
            return
        # Split the smiles string into a list and assign to the step
        mols[step] = smiles.split()
        bond_df = molecules_to_bond_energy_df(mols)
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
    entry = ttk.Entry(entry_frame, width=30)  # Adjust the width of the entry
    entry.grid(row=step_number + 3, column=0)
    steps_entries.append(entry)


# Function to add a random step with 1-3 small chemicals
# NOTE: WIP
# TODO: develop function
def generate_random_step():

    random_chemicals = []  # Extend list as needed
    random_step = ""
    chemletters = ["C", "N", "O", "F", "P", "S", "I", "c1cccc1", "c1ccccc1", "c1ccncc1"]
    for i in range(random.randint(1, 3)):
        random_step += "C" + random.choice(chemletters) + " "

    logger.info(f"Random step: {random_step}")

    step_number = len(steps_entries) + 1
    entry = ttk.Entry(entry_frame, width=30)  # Adjust the width of the entry
    entry.insert(0, f"{random_step}")
    entry.grid(row=step_number + 3, column=0)
    steps_entries.append(entry)


# Function to clear all existing steps and RDKit images
def reset():
    # use command line to close window and rerun the program 'python -m app.main'
    for widget in root.winfo_children():
        root.destroy()  # This will remove each widget
        os.system("python -m app.main")  # This will rerun the program


def on_scale_select(event):
    scale = scale_var.get()
    with open("util/dropdown_value.txt", "w") as f:
        f.write(scale)
        f.close()


# Read the SMILES guide from the text file
with open("smiles_notation_guide.txt", "r") as file:
    smiles_guide = file.read()

root = ThemedTk(theme="equilux")
root.title("Chemsong")
root.geometry("1200x900")  # Set a larger size for the window

# Title and Image at the top
top_frame = ttk.Frame(root)
top_frame.pack(fill=tk.X)

title_label = ttk.Label(top_frame, text="CHEMSONG", font=("Arial", 24))
title_label.pack(side=tk.LEFT, padx=10)

# Load and place the image
image = Image.open("img/chemsong.png")
image = image.resize((300, 300))  # Resize the image as needed
photo = ImageTk.PhotoImage(image)
image_label = tk.Label(top_frame, image=photo)
image_label.pack(side=tk.RIGHT, padx=10)

# Frame for entries and buttons
entry_frame = ttk.Frame(root)
entry_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.Y)

# Descriptive message
message_label = ttk.Label(
    entry_frame,
    text="Enter SMILES notation for each step, separate moelucules with a space:",
)
message_label.grid(row=0, column=0, padx=10, pady=5)


steps_entries = []
ttk.Button(entry_frame, text="Add Step", command=add_step).grid(
    row=1, column=0, padx=10, pady=5
)
ttk.Button(
    entry_frame, text="Process", command=lambda: process_reaction(steps_entries)
).grid(row=2, column=0, padx=10, pady=5)

random_button_text = "Generate Random Step"
ttk.Button(
    entry_frame, text=f"{random_button_text}", command=generate_random_step
).grid(row=3, column=0, padx=10, pady=5)

# Frame for SMILES guide
guide_frame = ttk.Frame(root)
guide_frame.pack(side=tk.RIGHT, padx=10, pady=10, fill=tk.Y)
guide_text = scrolledtext.ScrolledText(
    guide_frame, wrap=tk.WORD, width=40, height=20
)  # Adjust size
guide_text.insert(tk.INSERT, smiles_guide)
guide_text.config(state="disabled")
guide_text.pack()

# Play music button under the SMILES guide
play_button = ttk.Button(
    guide_frame, text="Play Music", command=lambda: play_notes(steps_entries)
)
play_button.pack(pady=10)
reset_button = ttk.Button(guide_frame, text="Reset", command=reset)
reset_button.pack(pady=10)

# Frame for RDKit images
image_frame = ttk.Frame(root)
image_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True, padx=10, pady=10)


# Creating the scale dropdown menu
ttk.Label(top_frame, text="Select musical scale:").pack(
    pady=5
)  # Adjust placement as needed
scale_var = tk.StringVar()
scale_menu = ttk.Combobox(top_frame, textvariable=scale_var, values=scale_list)
scale_menu.bind("<<ComboboxSelected>>", on_scale_select)
scale_menu.current("0")
scale_menu.pack(pady=5)  # Adjust placement as needed

root.mainloop()
