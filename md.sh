#!/usr/bin/env bash
set -e

#activate gromacs environment
# mamba activate gromacs

#path of the analysis
BASE="/mnt/wd_dnd/audrey_BT23/test"
WD="${BASE}/workdir"

#location and file name

    ##protein
P_PDBQT="${BASE}/pdbqt/protein/5WYZ.pdbqt"
PROTEIN="${WD}/00.protein"
INPUT_P="${PROTEIN}/protein.pdb"

    ##Ligand
PDBQT="${BASE}/pdbqt/ligand"
LIGAND="${WD}/00.ligand"
INPUT_L="${LIGAND}/"

    ##extra scripts
SCRIPTS="${BASE}/scripts"
export GMXLIB="${BASE}"

mkdir -p "${WD}" "${PROTEIN}" "${LIGAND}"

###############################################################
#----------- 0. Preping Input Files From Complex--------------#
###############################################################
#change from pdbqt to pdb

echo "Changing pdbqt format to pdb for protein..."

printf "open ${P_PDBQT}\nwrite format pdb #0 ${INPUT_P}\n" | chimera --nogui

echo "Splitting protein at chain gaps..."
python3 - "${INPUT_P}" << 'EOF'
import math, sys
pdb = sys.argv[1]
lines = open(pdb).readlines()
C, N = {}, {}
for l in lines:
    if l.startswith("ATOM"):
        key = (l[21], l[22:27])
        xyz = tuple(float(l[i:i+8]) for i in (30, 38, 46))
        name = l[12:16].strip()
        if name == "C": C[key] = xyz
        if name == "N": N[key] = xyz
order = list(dict.fromkeys((l[21], l[22:27]) for l in lines if l.startswith("ATOM")))
breaks = {b for a, b in zip(order, order[1:])
          if a[0] == b[0] and a in C and b in N and math.dist(C[a], N[b]) > 2.0}
out, seen = [], set()
for l in lines:
    if l.startswith("ATOM"):
        key = (l[21], l[22:27])
        if key in breaks and key not in seen:
            out.append("TER\n")
            print(f"  -> TER inserted before chain {key[0]} res {key[1].strip()}")
        seen.add(key)
    out.append(l)
open(pdb, "w").writelines(out)
EOF

echo "Changing pdbqt format to mol2 for ligand..."

for file in "${PDBQT}"/*.pdbqt; do
    [ -e "$file" ] || continue
    
    raw_name=$(basename "$file" .pdbqt)
    
    # 1. Remove '_minimized_out' from the back
    no_suffix="${raw_name%_minimized_out}"
    
    # 2. Remove everything up to the first underscore from the front (removes '5WYZ_')
    clean_name="${no_suffix#*_}"

    output_file="${INPUT_L}/${clean_name}.mol2"

    echo "Processing: $raw_name -> $clean_name.mol2"

    echo "Adding hydrogens with Chimera: $raw_name -> ${clean_name}.mol2"

    # Create temporary Chimera command file
    cat << EOF > prep_ligand.cmd
open ${file}
addh
write format mol2 #0 ${output_file}
stop
EOF

    # Run Chimera headless
    chimera --nogui prep_ligand.cmd > /dev/null 2>&1
    rm -f prep_ligand.cmd

    if [ -s "$output_file" ]; then
        echo "  -> Successfully saved ${output_file} with hydrogens."
    else
        echo "  -> ERROR: Chimera conversion failed for ${file}"
    fi
done

cd "${LIGAND}"
for ligand_file in *.mol2; do
    ###############################################################
    #-------0. Preparing Topology File for Protein & Ligand-------#
    ###############################################################
    # Check if any mol2 files exist to avoid running on an empty glob string
    [ -e "$ligand_file" ] || continue
    
    # Extract the ligand name without the .mol2 extension
    lig_name=$(basename "$ligand_file" .mol2)
    
    echo "Processing system: ${lig_name}"
    
    # Define system folders dynamically
    SYSTEM_DIR="${WD}/${lig_name}" 
    TOPOL_DIR="${SYSTEM_DIR}/01.Topology"
    SOLVATE_DIR="${SYSTEM_DIR}/02.solvate"
    
    # Create the subdirectories cleanly
    mkdir -p "${TOPOL_DIR}"
    mkdir -p "${SOLVATE_DIR}"

    echo "  -> Generating protein topology..."
    
    # Jump into the specific system topology folder to generate files cleanly
    cd "${TOPOL_DIR}"
    
    # FIX: Cleaned -ff flag to use name only, and localized outputs to current directory (.)
    gmx pdb2gmx \
        -f "${INPUT_P}" \
        -o "./protein_clean.gro" \
        -p "./topol.top" \
        -i "./posre.itp" \
        -ff "amber14sb_parmbsc1_cufix" \
        -water tip3p \
        -ignh
        
    echo "  -> System ${lig_name} topology folder complete!"
    echo "--------------------------------------------------------"
    
    # Move back to main ligand staging directory for the next iteration step
    cd "${LIGAND}"

    ###############################################################
    #---------1. Adding Hydrogen Atom and Ligand topology---------#
    ###############################################################
    
    # 1. Extract the actual molecule name from the 8th column of the first atom
    mol_name=$(awk '/@<TRIPOS>ATOM/ { getline; print $8; exit }' "$ligand_file")
    
    # 2. Replace the line immediately following @<TRIPOS>MOLECULE with the extracted name
    if [ -n "$mol_name" ]; then
        sed -i '/@<TRIPOS>MOLECULE/{n;s/.*/'"${mol_name}"'/}' "$ligand_file"
        echo "  -> Updated molecule name in ${ligand_file} to ${mol_name}"
    else
        echo "  -> Warning: Could not extract molecule name for ${ligand_file}"
    fi
    echo "Runnning sort_mol2_bonds.pl from Dr. Lemkul"
    perl ${SCRIPTS}/sort_mol2_bonds.pl "${ligand_file}" "${SYSTEM_DIR}/${lig_name}_fix.mol2"
    echo "Script done, moving on to next conversion tool to itp format!"
    # move to each ligands directory to run acpype
    (cd "${SYSTEM_DIR}" && mamba run -n acype acpype -i "${lig_name}_fix.mol2" -a gaff2 -c bcc)

    #edit config of the file
    gmx editconf -f "${SYSTEM_DIR}/${lig_name}_fix.acpype/${lig_name}_fix_GMX.gro" -o "${TOPOL_DIR}/unl.gro"

    #adding changes to top.itp
    # Copy the GMX itp file into 01.Topology
    # 1. Copy ACPYPE's native ITP file
    cp "${SYSTEM_DIR}/${lig_name}_fix.acpype/${lig_name}_fix_GMX.itp" "${TOPOL_DIR}/${lig_name}_fix_GMX.itp"

    # 2. Extract [ atomtypes ] into its own file
    awk '
        /\[[[:space:]]*atomtypes[[:space:]]*\]/ { show=1; print; next }
        show && /^\[/ { show=0 }
        show { print }
    ' "${TOPOL_DIR}/${lig_name}_fix_GMX.itp" > "${TOPOL_DIR}/${lig_name}_atomtypes.itp"

    # 3. Strip [ atomtypes ], [ defaults ], [ system ], and [ molecules ] out of the main ITP
    awk '
        /\[[[:space:]]*(atomtypes|defaults|system|molecules)[[:space:]]*\]/ { skip=1; next }
        skip && /^\[/ { skip=0 }
        !skip { print }
    ' "${TOPOL_DIR}/${lig_name}_fix_GMX.itp" > "${TOPOL_DIR}/temp.itp" && mv "${TOPOL_DIR}/temp.itp" "${TOPOL_DIR}/${lig_name}_fix_GMX.itp"

    # 4. Extract the moleculetype name
    mol_name=$(awk '/\[ moleculetype \]/ {
        while (getline > 0) {
            if ($0 ~ /^;/ || $0 ~ /^[[:space:]]*$/) continue;
            print $1;
            exit;
        }
    }' "${TOPOL_DIR}/${lig_name}_fix_GMX.itp")

    # 5. Insert atomtypes inclusion right after forcefield.itp
    sed -i '/forcefield\.itp"/a\; Include ligand atomtypes\
#include "'"${TOPOL_DIR}/${lig_name}_atomtypes.itp"'"' "${TOPOL_DIR}/topol.top"

    # 6. Insert ligand topology inclusion right before [ system ]
        sed -i '/\[[[:space:]]*system[[:space:]]*\]/i\; Include ligand topology\
#include "'"${TOPOL_DIR}/${lig_name}_fix_GMX.itp"'"\
' "${TOPOL_DIR}/topol.top"

    # 7. Append ligand name to the [ molecules ] list
    echo "${mol_name}             1" >> "${TOPOL_DIR}/topol.top"

    echo "  -> Successfully linked ${lig_name} topology into topol.top!"
    echo "--------------------------------------------------------"
    
    cd "${TOPOL_DIR}"
    TOPOL_TOP="${TOPOL_DIR}/topol.top"
    LIGAND_GRO="${SYSTEM_DIR}/${lig_name}_fix.acpype/${lig_name}_fix_GMX.gro"
    COMPLEX_GRO="${TOPOL_DIR}/complex.gro"
    prot_atoms=$(sed -n '2p' protein_clean.gro | tr -d ' ')
    lig_atoms=$(sed -n '2p' "${LIGAND_GRO}" | tr -d ' ')
    total_atoms=$((prot_atoms + lig_atoms))

    # 2. Write Header (Line 1)
    echo "Complex of protein and ${lig_name}" > "${COMPLEX_GRO}"

    # 3. Write New Total Atom Count (Line 2)
    echo " ${total_atoms}" >> "${COMPLEX_GRO}"

    # 4. Append Protein Atom Lines (skip first 2 lines & last line)
    sed '1,2d;$d' protein_clean.gro >> "${COMPLEX_GRO}"
    echo "Finish append protein atom lines"

    # 5. Append Ligand Atom Lines (skip first 2 lines & last line)
    sed '1,2d;$d' "${LIGAND_GRO}" >> "${COMPLEX_GRO}"
    echo "Finish ligand atom lines"

    # 6. Append Box Vectors (take last line from protein gro)
    tail -n 1 protein_clean.gro >> "${COMPLEX_GRO}"

    echo "  -> Created ${COMPLEX_GRO} with ${total_atoms} total atoms!"

    gmx editconf -f complex.gro -o "${SOLVATE_DIR}"/newbox.gro -bt dodecahedron -d 1.0

    ###############################################################
    #---------------------2. Solvate Complex----------------------#
    ###############################################################
    ION_DIR="${SYSTEM_DIR}/03.ion"
    
    mkdir -p "${ION_DIR}"

    echo " doing solvation"
    gmx solvate -cp "${SOLVATE_DIR}/newbox.gro" -cs spc216.gro -p "${TOPOL_TOP}" -o "${SOLVATE_DIR}/solv.gro"

    cd "${ION_DIR}"
    echo "setting up ion addition"
    gmx grompp -f "${SCRIPTS}/ions.mdp" -c "${SOLVATE_DIR}/solv.gro" -p "${TOPOL_TOP}" -o ./ions.tpr
    echo "adding ion"
    echo "SOL" | gmx genion -s ./ions.tpr -o "${SOLVATE_DIR}/solv_ions.gro" -p "${TOPOL_TOP}" -pname NA -nname CL -neutral

    ###############################################################
    #-------------------3. Energy Minimization--------------------#
    ###############################################################
    EM_DIR="${SYSTEM_DIR}/04.em"
    
    mkdir -p "${EM_DIR}"

    # Minimize the cleaned structure
    echo "doing minimization"
    cd "${EM_DIR}"
    gmx grompp -f "${SCRIPTS}/em.mdp" -c "${SOLVATE_DIR}/solv_ions.gro" -r "${SOLVATE_DIR}/solv_ions.gro" -p "${TOPOL_TOP}" -o ./em.tpr
    gmx mdrun -v -deffnm ./em
    echo "finish minimization"


    ###############################################################
    #----------------------4. Equilibration-----------------------#
    ###############################################################
    EQU_DIR="${SYSTEM_DIR}/05.equ"
    mkdir -p "${EQU_DIR}"
    cd "${EQU_DIR}"

    echo "restaining ligand"

    #Indexing for equilibration
    printf "0 & ! a H*\nq\n" | gmx make_ndx -f "${TOPOL_DIR}/unl.gro" -o ./index_unl1.ndx
    echo "3" | gmx genrestr -f "${TOPOL_DIR}/unl.gro" -n ./index_unl1.ndx -o ./posre_unl1.itp -fc 1000 1000 1000
    
    # Ensure posre_unl.itp is in the topology directory
    cp "${EQU_DIR}/posre_unl1.itp" "${TOPOL_DIR}/posre_unl1.itp" 2>/dev/null || true

    # Inject POSRES_LIG block into topol.top above water topology
    python3 -c '
topol_path = "'"${TOPOL_TOP}"'"
with open(topol_path) as f:
    lines = f.readlines()
if not any("POSRES_LIG" in l for l in lines):
    out = []
    for line in lines:
        if line.strip().startswith("[ system ]"):
            out += ["; Ligand position restraints\n", "#ifdef POSRES_LIG\n",
                    "#include \"posre_unl1.itp\"\n", "#endif\n\n"]
        out.append(line)
    with open(topol_path, "w") as f:
        f.writelines(out)
'

    #Equilibration
    cd "${TOPOL_DIR}"
    echo "equilibrating NVT"
    gmx grompp -f "${SCRIPTS}/nvt.mdp" -c "${EM_DIR}/em.gro" -r "${EM_DIR}/em.gro" -p ./topol.top -o "${EQU_DIR}/nvt.tpr"
    
    echo "equilibration of NVT done, now running the simulation"
    gmx mdrun -deffnm "${EQU_DIR}/nvt" -update cpu -ntomp 8
    echo "NVT simulation complete"
    echo "Temperature" | gmx energy -f "${EQU_DIR}/nvt.edr" -o "${EQU_DIR}/temperature.xvg"
    
    echo "start NPT"
    gmx grompp -f "${SCRIPTS}/npt.mdp" -c "${EQU_DIR}/nvt.gro" -t "${EQU_DIR}/nvt.cpt" -r "${EQU_DIR}/nvt.gro" -p ./topol.top -o "${EQU_DIR}/npt.tpr"
    gmx mdrun -deffnm "${EQU_DIR}/npt" -update cpu -pin on -ntomp 8
    echo "Pressure"    | gmx energy -f "${EQU_DIR}/npt.edr" -o "${EQU_DIR}/pressure.xvg"
    echo "Density"     | gmx energy -f "${EQU_DIR}/npt.edr" -o "${EQU_DIR}/density.xvg"
    echo "NPT done"

    ###############################################################
    #-------------------------5. Run MDS--------------------------#
    ###############################################################
    MD_DIR="${SYSTEM_DIR}/06.md"
    mkdir -p "${MD_DIR}"
    cd "${MD_DIR}"

    echo "Run MD!"
    gmx grompp -f "${SCRIPTS}/md.mdp" -c "${EQU_DIR}/npt.gro" -t "${EQU_DIR}/npt.cpt" -p "${TOPOL_TOP}" -o ./md.tpr
    gmx mdrun -deffnm ./md -update cpu -ntomp 8
    cd "${LIGAND}"

done

echo "Script complete! Directory tree and topologies generated successfully."