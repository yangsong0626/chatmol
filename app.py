from flask import Flask, render_template, request, jsonify, send_file
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Crippen
from rdkit.Chem.rdchem import Mol
import base64
from io import BytesIO, StringIO
import random

app = Flask(__name__)

# List of functional groups to add
FUNCTIONAL_GROUPS = [
    "C(=O)O",  # Carboxylic acid
    "C(=O)N",  # Amide
    "C(=O)C",  # Ketone
    "NC",      # Methylated amine
    "O",       # Hydroxyl
    "OC",      # Methoxy
    "Cl",      # Chlorine
    "Br",      # Bromine
    "F",       # Fluorine
    "C(F)(F)F",#CF3
    "C#N",     # Nitrile
    "C1CCCCC1",# Cyclohexane
    "C1CC1",   # Cyclopropyl
]

# Function to generate analogs by joining functional groups
def generate_analogs(smiles, num_analogs=100):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    analogs = []
    for _ in range(num_analogs):
        # Randomly select a functional group
        fg = random.choice(FUNCTIONAL_GROUPS)
        fg_mol = Chem.MolFromSmiles(fg)
        if not fg_mol:
            continue

        # Create a writable copy of the input molecule
        new_mol = Chem.RWMol(mol)

        # Combine the input molecule and the functional group
        combined_mol = Chem.CombineMols(new_mol, fg_mol)

        # Select random atoms to form a bond
        mol_atom_idx = random.randint(0, new_mol.GetNumAtoms() - 1)  # Atom in the input molecule
        fg_atom_idx = random.randint(0, fg_mol.GetNumAtoms() - 1)    # Atom in the functional group

        # Create an editable molecule
        editable_mol = Chem.EditableMol(combined_mol)

        # Add a bond between the selected atoms
        editable_mol.AddBond(
            mol_atom_idx,  # Atom in the input molecule
            new_mol.GetNumAtoms() + fg_atom_idx,  # Atom in the functional group
            order=Chem.rdchem.BondType.SINGLE  # Bond type (single bond)
        )

        # Convert the editable molecule back to a regular molecule
        new_analog = editable_mol.GetMol()

        # Sanitize the molecule to ensure it's valid
        try:
            Chem.SanitizeMol(new_analog)
            analogs.append(new_analog)
        except:
            continue  # Skip invalid molecules

    return analogs

# Function to calculate logD
def calculate_logD(mol):
    return Crippen.MolLogP(mol)

# Function to convert molecule to image (base64)
def mol_to_image_base64(mol):
    img = Draw.MolToImage(mol)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode("utf-8")

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/generate', methods=['POST'])
def generate():
    smiles = request.form['smiles']
    analogs = generate_analogs(smiles)
    if not analogs:
        return jsonify({"error": "Invalid SMILES string"}), 400

    # Convert molecules to images and calculate logD
    input_image = mol_to_image_base64(Chem.MolFromSmiles(smiles))
    analog_data = []
    for mol in analogs:
        analog_data.append({
            "image": mol_to_image_base64(mol),
            "smiles": Chem.MolToSmiles(mol),
            "logD": calculate_logD(mol)
        })

    return jsonify({
        "input_image": input_image,
        "analogs": analog_data
    })

@app.route('/download', methods=['POST'])
def download():
    selected_smiles = request.form.getlist('smiles[]')  # Get selected SMILES strings
    if not selected_smiles:
        return jsonify({"error": "No analogs selected"}), 400

    # Create an SDF file
    sdf_data = StringIO()
    writer = Chem.SDWriter(sdf_data)
    for smiles in selected_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            writer.write(mol)
    writer.close()

    # Save the SDF file temporarily
    sdf_filename = "selected_analogs.sdf"
    with open(sdf_filename, "w") as f:
        f.write(sdf_data.getvalue())

    return send_file(sdf_filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)