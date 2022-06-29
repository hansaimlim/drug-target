"""Update templates for the given pdb id."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from pathlib import Path
from subprocess import Popen
import os
import sys
import argparse
import logging
import gzip
import json
import numpy as np
from django.conf import settings
import biotite
from biotite.structure.io import pdbx
from biotite.structure.io import pdb
from biotite import structure as Structure
from . import rosetta_dock as rosetta
from . import protein
from . import utils
from .utils import download_pdb, run_eppic, str2bool
from .utils import extract_chain_info_from_cif, save_json, compute_euclidean_distance


class Residue:
    """Represent each amino acid residue from biotite.
    
    This class is devised to solve issues with intertion codes in cif files.
    self.id is an internally used sequential integer id, different from res_id
    self.id is unique per residue, matching to the residue_starts obtained from biotite.structure.get_residue_starts
    """

    def __init__(self,
        id=None, atomarray=None, sasa=None, chain_id=None,
        relative_sasa=None, res_name=None, res_id=None
    ):
        self.id = id # internal ID. sequential integer in a chain; may differ from res_id
        self.atomarray = atomarray
        self.sasa = sasa
        self.chain_id = chain_id
        self.res_id = res_id
        self.ins_code = None
        self.res_name = res_name
        self.CA = None
        self.com = None # center of mass
        self.pair_potential = 0
        self.relative_sasa = relative_sasa
        self.is_hotspot = False # default False, updated later to True if conditions met

        if self.atomarray is None:
            raise ValueError("Atom array for a residue cannot be None.")
        if self.id is None:
            raise ValueError("Residue internal ID cannot be None")

        if self.chain_id is None:
            self.chain_id = self.atomarray[0].chain_id
        if self.CA is None:
            self.CA = self.atomarray[self.atomarray.atom_name == 'CA']
        if self.res_id is None:
            self.res_id = self.atomarray[0].res_id
        if self.res_name is None:
            self.res_name = self.atomarray[0].res_name

    def center_of_mass(self):
        """Return center of mass (com) for this residue.
        If com is None, calculate and return the com.
        """

        if self.com is None:
            self.com = protein.compute_center_of_mass_residue(self.atomarray, res_name=self.res_name)
        return self.com


class Residues:
    """Define a set of Residue objects.
    
    get_residue_by_id
    get_all_residue_ids
    get_hotspot_ids
    atom_count
    """

    def __init__(self,
        atomarray, sasaarray
    ):
        self.atomarray = atomarray
        self.sasaarray = sasaarray
        self.residue_starts = Structure.get_residue_starts(atomarray, add_exclusive_stop=True)
        self.residues = {}
        self.last_residue_id = 0
        stasa = protein.AminoAcid['standardASA']
    
        for i in range(len(self.residue_starts[:-1])):
            rs = self.residue_starts[i] # residue start
            re = self.residue_starts[i+1] # residue end
            atoms = atomarray[rs:re]
            sasa = self.sasaarray[rs]
            res_name = atoms[0].res_name
            if (res_name not in stasa) or (np.isnan(sasa)):
                rel_sasa = 999
            else:
                rel_sasa = (100.0 * sasa) / stasa[res_name]
            residue = Residue(id=self.last_residue_id+1, 
                chain_id=atoms[0].chain_id,
                atomarray=atoms, sasa=sasa, relative_sasa=rel_sasa,
                res_name=res_name, res_id=atoms[0].res_id
                )
            self.residues[self.last_residue_id+1] = residue
            self.last_residue_id += 1
    
    def __len__(self):
        return len(self.residues)

    def residue_count(self):
        return len(self.residues)

    def atom_count(self):
        if self.atomarray is None:
            return 0
        else:
            return len(self.atomarray)
    
    def get_residue_by_id(self, id):
        if id in self.residues:
            return self.residues[id]
        else:
            return None
    
    def get_all_residue_ids(self):
        if len(self.residues) > 0:
            return sorted(list(self.residues.keys()))
        else:
            return []

    def get_hotspot_ids(self):
        hotspots = sorted([self.residues[k].id for k in self.get_all_residue_ids() if self.residues[k].is_hotspot])
        return hotspots


def residues_in_interaction_chain(residues1, residues2):
    """Compute interacting residues from two chains.
    
    Returns two lists, containing interface residue ids, one for each chain.
    residues1 and residues2 are Residues objects for each chain. Each residue will be processed by
    residues_in_interaction method above.
    """
    residue_ids_1 = residues1.get_all_residue_ids()
    residue_ids_2 = residues2.get_all_residue_ids()
    interacting_residues = {
        1: [],
        2: []
    }
    for id1 in residue_ids_1:
        residue1 = residues1.get_residue_by_id(id1)
        atoms1 = residue1.atomarray
        if id1 in interacting_residues[1]:
            continue

        for id2 in residue_ids_2:
            if id2 in interacting_residues[2]:
                continue
            residue2 = residues2.get_residue_by_id(id2)
            atoms2 = residue2.atomarray

            if protein.residues_in_interaction(atoms1, atoms2):
                interacting_residues[1].append(id1)
                interacting_residues[2].append(id2)
                break
    for k in interacting_residues:
        interacting_residues[k] = sorted(list(set(interacting_residues[k])))
    
    return interacting_residues


def get_nearby_residues(residues, interacting_residue_ids, threshold=6.0):
    """Collect nearby residues of interacting residues on a single chain."""

    residue_ids = residues.get_all_residue_ids()
    nearby_residue_ids = []
    for id1 in interacting_residue_ids:
        residue1 = residues.get_residue_by_id(id1)
        atoms1 = residue1.atomarray

        for id2 in residue_ids:
            if id1 == id2:
                continue
            if id2 in nearby_residue_ids:
                continue
            if id2 in interacting_residue_ids:
                continue
            residue2 = residues.get_residue_by_id(id2)
            atoms2 = residue2.atomarray
            if protein.residues_in_contact(atoms1, atoms2, distance_threshold=threshold, return_distance=False):
                nearby_residue_ids.append(id2)

    return nearby_residue_ids

def compute_residue_pair_potential(
        residues, interface_residue_ids, residues2, com_dist_max=7.0, res_id_diff=4,
        min_potential=18.0, max_sasa=20.0
        ):
    """Compute residue pair potential for each residue in interface.
    interface_residue_ids: list of residue ids in interface; [interacting_residue_ids] + [nearby_residue_ids]
    For each pair of residues:
        if |res_id1 - res_id2| >= res_id_diff  # sequentially distant
        and if com_dist(res1, res2) <= com_dist_max  # spatially close
        compute pairpot(res1, res2)
    Returns list of (res_id, 3-letter-code, sum-pair-potential)
    """
    residue_ids = residues.get_all_residue_ids()
    residue_ids2 = residues2.get_all_residue_ids()
    residue_pair_potentials = []
    for id1 in interface_residue_ids:
        residue1 = residues.get_residue_by_id(id1)
        chain_id = residue1.chain_id
        res_name = residue1.res_name
        com1 = residue1.center_of_mass()
        resno1 = residue1.res_id  # residue number in chain
        sasa = residue1.sasa
        rel_sasa = residue1.relative_sasa
        total_pot = 0.0

        for id2 in residue_ids:
            residue2 = residues.get_residue_by_id(id2)
            resno2 = residue2.res_id
            if np.abs(resno1 - resno2) >= res_id_diff:
                # sequentially distant
                com2 = residue2.center_of_mass()

                com_dist = protein.compute_euclidean_distance(com1, com2)
                if com_dist <= com_dist_max:
                    #spatially close
                    pot = protein.ResiduePairPotential(res_name, residue2.res_name)
                    total_pot += pot
        for id2 in residue_ids2:
            residue2 = residues2.get_residue_by_id(id2)
            com2 = residue2.center_of_mass()

            com_dist = protein.compute_euclidean_distance(com1, com2)
            if com_dist <= com_dist_max:
                # spatially close
                pot = protein.ResiduePairPotential(res_name, residue2.res_name)
                total_pot += pot
        
        residue1.pair_potential = total_pot
        if np.abs(residue1.pair_potential) >= min_potential:
            if rel_sasa <= max_sasa:
                residue1.is_hotspot = True

        residue_pair_potentials.append({
            'id': id1,
            'chain_id': chain_id,
            'res_name': res_name,
            'total_potential': total_pot,
            'relative_sasa': rel_sasa,
            'sasa': sasa,
            'res_id': resno1,
            'ins_code': residue1.ins_code,
            'is_hotspot': residue1.is_hotspot
        })
    
    return residue_pair_potentials


def collect_hotregion(residues1, hotspots1, residues2, hotspots2, threshold=6.5):
    """Collect a network of hot spots, named hot region.
    
    hotspots (1 or 2) are list of residue internal ids specific to the residues (1 or 2)
    hotregion is a dict, (hotspot) res_id -> list of other hotspot res_ids within hotregion threshold
    endo1: pairs of hot spots within chain 1
    endo2: pairs of hot spots within chain 2
    exo: pairs of hot spots across chains
    """

    hotregions = {
        'endo1': [],
        'endo2': [],
        'exo': []
    }
    for id1 in hotspots1:
        residue1 = residues1.get_residue_by_id(id1)
        com1 = residue1.center_of_mass()

        for id1_2 in hotspots1:
            # collect endo0
            residue1_2 = residues1.get_residue_by_id(id1_2)
            if id1 == id1_2:
                continue
            com1_2 = residue1_2.center_of_mass()
            dist = compute_euclidean_distance(com1, com1_2)
            if dist <= threshold:
                hotregions['endo1'].append((id1, id1_2))

        for id2 in hotspots2:
            # collect exo
            residue2 = residues2.get_residue_by_id(id2)
            com2 = residue2.center_of_mass()
            dist = compute_euclidean_distance(com1, com2)
            if dist <= threshold:
                hotregions['exo'].append((id1, id2))
    
    for id2 in hotspots2:
        residue2 = residues2.get_residue_by_id(id2)
        com2 = residue2.center_of_mass()

        for id2_2 in hotspots2:
            #collect endo1
            if id2 == id2_2:
                continue
            residue2_2 = residues2.get_residue_by_id(id2_2)
            com2_2 = residue2_2.center_of_mass()
            dist = compute_euclidean_distance(com2, com2_2)
            if dist <= threshold:
                hotregions['endo2'].append((id2, id2_2))
    return hotregions

def collect_interfaces_from_cif(
    filename, assembly_id=None, chain_mode='specific', chain1='', chain2='',
    sasa_probe_radius=1.4, sasa_point_number=1000, sasa_ignore_ions=True,
    scaffold_distance_threshold=6.0,
    hotspot_com_distance_threshold=7.0,
    hotspot_res_id_diff=4,
    hotspot_min_total_potential=18.0,
    hotspot_max_sasa=20.0,
    hotregion_com_distance_threshold=6.5,
    combined_file_output_dir=None,
    clean=True
    ):
    """Collect interfaces by reading cif file."""
    """Collect chain/atom info from one interface file using biotite.
    In case of .pdb file, set assembly_id=None. A .pdb file with assembly_id will return an exception
    Set chain_mode="dataprep" to process pointcloud dataset and skip unnecessary parts
    """
    if assembly_id is None:
        # .pdb file goes here
        cifStructure = protein.load_cif(filename, download=False, structure_type="structure")
    else:
        # .cif file only
        cifStructure = protein.load_assembly(filename, assembly_id)
    is_valid = False

    try:
        chain_ids = Structure.get_chains(cifStructure)
    except AttributeError:
        error_message = "ATTRERROR: {} Could not be read properly. Chains {} - {}.".format(
            filename, chain1, chain2
        )
        return (False, error_message)


    if chain_mode.lower() == 'specific':
        if chain1 not in chain_ids:
            error_message = "Chain ID ({}) not found in the input file.".format(chain1)
            if clean:
                Popen("rm -rf {}".format(filename), shell=True).wait()
            return (is_valid, error_message)
        elif chain2 not in chain_ids:
            error_message = "Chain ID ({}) not found in the input file.".format(chain2)
            if clean:
                Popen("rm -rf {}".format(filename), shell=True).wait()
            return (is_valid, error_message)
        else:
            chain_pairs = [(chain1, chain2)]
    
    # save a file for the two chains
    atomarray1 = cifStructure[cifStructure.chain_id == chain1]
    atomarray2 = cifStructure[cifStructure.chain_id == chain2]
    combined_atomarray = atomarray1 + atomarray2
    combinedFileObj = pdbx.PDBxFile()
    pdbx.set_structure(combinedFileObj, combined_atomarray, data_block="structure")
    combinedFileObj.set_category("citation", {"title":" "})
    if combined_file_output_dir is None:
        combined_file_output_dir = settings.HMI_PATHS['TEMPORARY_STRUCTURE_PATH']

    combined_file = os.path.join(combined_file_output_dir, utils.random_filename(size=14))
    combinedFileObj.write(combined_file)

    interface = {
        'pdb_file': combined_file,
    }

    # compute sasa for each chain
    # compute c.o.m for each chain
    # collect interface (within vdw1+vdw2+0.5 and scaffold Ca within 6.0)
    # compute pair-potential for each residue in interface
    # collect hot spots
    #   - c.o.m. dist <= 7.0 and dResNo >= 4
    #   - relCompASA <= 20.0% and abs(sum(pairpotential)) >= 18.0 => hot spot
    # compute Ca distance for each hot spot pair
    #   - hot region network - residues within Ca dist 6.5 
    Chains = {}
    chain_starts = Structure.get_chain_starts(combined_atomarray, add_exclusive_stop=True)

    if chain_mode.lower() == 'dataprep':
        sasa_per_atom, sasa_per_res = protein.compute_sasa(combined_atomarray,
            probe_radius=sasa_probe_radius, point_number=sasa_point_number, ignore_ions=sasa_ignore_ions)
        sasa1 = sasa_per_res[:chain_starts[1]]
        sasa2 = sasa_per_res[chain_starts[1]:]
        residues1 = Residues(atomarray1, sasa1)
        residues2 = Residues(atomarray2, sasa2)
    else:

        try:
            sasa_per_atom, sasa_per_res = protein.compute_sasa(combined_atomarray,
                probe_radius=sasa_probe_radius, point_number=sasa_point_number, ignore_ions=sasa_ignore_ions)
            sasa1 = sasa_per_res[:chain_starts[1]]
            sasa2 = sasa_per_res[chain_starts[1]:]
            residues1 = Residues(atomarray1, sasa1)
            residues2 = Residues(atomarray2, sasa2)

        except:
            error_message = "Structures could not be parsed properly. Please check the inpuf file format."
            if clean:
                Popen("rm -rf {} {}".format(filename, combined_file), shell=True).wait()
            return (False, error_message)

    interacting_residues = residues_in_interaction_chain(residues1, residues2)
    nearby_residues1 = get_nearby_residues(residues1, interacting_residues[1], threshold=scaffold_distance_threshold)
    nearby_residues2 = get_nearby_residues(residues2, interacting_residues[2], threshold=scaffold_distance_threshold)
    interface_residues1 = sorted(list(set(interacting_residues[1] + nearby_residues1)))
    interface_residues2 = sorted(list(set(interacting_residues[2] + nearby_residues2)))

    residue_pairpot1 = compute_residue_pair_potential(
        residues1, interface_residues1, residues2, com_dist_max=hotspot_com_distance_threshold, res_id_diff=hotspot_res_id_diff,
        min_potential=hotspot_min_total_potential, max_sasa=hotspot_max_sasa
        )
    residue_pairpot2 = compute_residue_pair_potential(
        residues2, interface_residues2, residues1, com_dist_max=hotspot_com_distance_threshold, res_id_diff=hotspot_res_id_diff,
        min_potential=hotspot_min_total_potential, max_sasa=hotspot_max_sasa
        )
    hotspots1 = residues1.get_hotspot_ids()
    hotspots2 = residues2.get_hotspot_ids()
    hotregions = collect_hotregion(residues1, hotspots1, residues2, hotspots2, threshold=hotregion_com_distance_threshold)

    interface['is_valid'] = True
    interface['chain1'] = {
        'chain_id': chain1,
        'residues': residues1,
        'interacting_residue_ids': interacting_residues[1],
        'nearby_residue_ids': nearby_residues1,
        'interface_residue_ids': interface_residues1,
        'residue_pair_potentials': residue_pairpot1,
        'hotspot_residue_ids': hotspots1,
    }
    interface['chain2'] = {
        'chain_id': chain2,
        'residues': residues2,
        'interacting_residue_ids': interacting_residues[2],
        'nearby_residue_ids': nearby_residues2,
        'interface_residue_ids': interface_residues2,
        'residue_pair_potentials': residue_pairpot2,
        'hotspot_residue_ids': hotspots2,
    }
    interface['hotregions'] = hotregions

    return (True, interface)


def check_downloaded_cif(filepath):
    """Open cif and check validity"""
    try:
        obj = pdbx.PDBxFile.read(filepath)
        structure = pdbx.get_assembly(obj, altloc="occupancy")
        if structure is None:
            is_valid = False
        else:
            is_valid = True
    except:
        is_valid = False
    return is_valid


def get_rosetta_scores(
        twochain_file, chain_ids, interfaceDir,
        rosetta_prepack=None, rosetta_db=None, rosetta_dock=None
    ):
    if rosetta_prepack is None:
        rosetta_prepack = settings.HMI_TOOLS['ROSETTAPREPACK']
    if rosetta_db is None:
        rosetta_db = settings.HMI_TOOLS['ROSETTA_DB']
    if rosetta_dock is None:
        rosetta_dock = settings.HMI_TOOLS['ROSETTADOCK']

    try:
        prepackedfile, scorefile_ = rosetta.run_rosettaprepack(
                twochain_file, chain_ids, interfaceDir,
                rosetta_prepack, rosetta_db)
        dockfile, scorefile = rosetta.run_rosettadock(
            prepackedfile, chain_ids, interfaceDir, 
            rosetta_dock, rosetta_db)
        
        t_score, i_score = rosetta.parse_rosettascore(scorefile)

        Popen("rm -f {} {} {} {}".format(
            prepackedfile, scorefile_, scorefile, dockfile
        ), shell=True).wait()
        rosetta_computed = "True"
    except:
        t_score = 0.0
        i_score = 0.0
        rosetta_computed = "False"
    output = {
        'interaction_score': i_score,
        'total_score': t_score,
        'rosetta_computed': rosetta_computed,
        }
    return output

def check_interface_files(interfaceDir, pdbid, chain1, chain2):
    twochain_file = os.path.join(interfaceDir, "{}_{}_{}.cif".format(pdbid, chain1, chain2))
    chain1_file = os.path.join(interfaceDir, "{}_{}.cif".format(pdbid, chain1))
    chain2_file = os.path.join(interfaceDir, "{}_{}.cif".format(pdbid, chain2))
    interface_residue_file = os.path.join(interfaceDir, "{}_{}_{}.intres".format(pdbid, chain1, chain2))
    interface_info_file = os.path.join(interfaceDir, "{}_{}_{}.json".format(pdbid, chain1, chain2))
    files = {
        'twochain_file': twochain_file,
        'chain1_file': chain1_file,
        'chain2_file': chain2_file,
        'interface_residue_file': interface_residue_file,
        'interface_info_file': interface_info_file,
    }
    for fk in files:
        fn = files[fk]
        if not os.path.exists(fn):
            return (False, fk)
    return (True, '')

def write_interface_files(interface, interfaceDir, pdbid, pair, chain_info, assembly_id=None, HMI_TOOLS=None, validate=True):
    if HMI_TOOLS is None:
        HMI_TOOLS = settings.HMI_TOOLS
    twochain_file = os.path.join(interfaceDir, "{}_{}_{}.cif".format(pdbid, pair[0], pair[1]))
    chain1_file = os.path.join(interfaceDir, "{}_{}.cif".format(pdbid, pair[0]))
    chain2_file = os.path.join(interfaceDir, "{}_{}.cif".format(pdbid, pair[1]))
    interface_residue_file = os.path.join(interfaceDir, "{}_{}_{}.intres".format(pdbid, pair[0], pair[1]))
    interface_info_file = os.path.join(interfaceDir, "{}_{}_{}.json".format(pdbid, pair[0], pair[1]))

    chain1_info = chain_info[pair[0]] # ('P20231', 'TRYB2_HUMAN', '9606', 'Homo sapiens', 'human')
    chain2_info = chain_info[pair[1]]

    residues1 = interface['chain1']['residues']
    chain_id1 = interface['chain1']['chain_id']
    interface_residues1 = interface['chain1']['interface_residue_ids']
    hotspot_residues1 = residues1.get_hotspot_ids()

    residues2 = interface['chain2']['residues']
    chain_id2 = interface['chain2']['chain_id']
    interface_residues2 = interface['chain2']['interface_residue_ids']
    hotspot_residues2 = residues2.get_hotspot_ids()

    # filtering
    # at least 1 hot spot and 15 interface residues for the mimicked chain
    interface_is_valid = False
    if chain1_info[2] == '9606' and chain2_info[2] == '9606':
        #endogeneous interface
        if len(interface_residues1) >= 15 and len(hotspot_residues1) > 0:
            interface_is_valid = True
        if len(interface_residues2) >= 15 and len(hotspot_residues2) > 0:
            interface_is_valid = True                
    elif chain1_info[2] == '9606':
        #chain 2 will be mimicked
        if len(interface_residues2) >= 15 and len(hotspot_residues2) > 0:
            interface_is_valid = True
    elif chain2_info[2] == '9606':
        #chain 1 will be mimicked
        if len(interface_residues1) >= 15 and len(hotspot_residues1) > 0:
            interface_is_valid = True
    if validate and (not interface_is_valid):
        # validation unnecessary for pointcloud dataprep
        return False

    pdb_file = interface['pdb_file']
    # move temporary two-chain file to permanent filename
    Popen("mv {} {}".format(pdb_file, twochain_file),shell=True).wait()
    logging.debug("Moved {} to {}".format(pdb_file, twochain_file))

    chain_ids = "{}_{}".format(pair[0],pair[1])
    rosetta_output = get_rosetta_scores(
                        twochain_file, chain_ids, interfaceDir,
                        rosetta_prepack=HMI_TOOLS['ROSETTAPREPACK'], 
                        rosetta_db=HMI_TOOLS['ROSETTA_DB'], 
                        rosetta_dock=HMI_TOOLS['ROSETTADOCK']
                    )    

    obj1 = pdbx.PDBxFile()
    pdbx.set_structure(obj1, residues1.atomarray, data_block="structure")
    obj1.set_category("citation", {"title": " "})
    obj1.write(chain1_file)
    obj2 = pdbx.PDBxFile()
    pdbx.set_structure(obj2, residues2.atomarray, data_block="structure")
    obj2.set_category("citation", {"title": " "})
    obj2.write(chain2_file)
    
    interfaceInfo = {
        "InterfaceName": "{}_{}_{}".format(pdbid, chain_id1, chain_id2),
        "ChainID_left": chain_id1,
        "ChainID_right": chain_id2,
        "PDBID": pdbid,
        "InterfaceLength_left": len(interface_residues1),
        "InterfaceLength_right": len(interface_residues2),
        "HotspotCount_left": len(hotspot_residues1),
        "HotspotCount_right": len(hotspot_residues2),
        "InterfaceFile_left": str(chain1_file),
        "InterfaceFile_right": str(chain2_file),
        "UniprotAccession_left": 'None' if chain1_info[0]==pdbid else chain1_info[0],
        "Genename_left": 'None' if chain1_info[1]==pdbid else chain1_info[1],
        "is_human_left": "True" if chain1_info[2] == '9606' else "False",
        "UniprotAccession_right": 'None' if chain2_info[0]==pdbid else chain2_info[0],
        "Genename_right": 'None' if chain2_info[1]==pdbid else chain2_info[1],
        "is_human_right": "True" if chain2_info[2] == '9606' else "False",
        "ChainLength_left": len(residues1),
        "ChainLength_right": len(residues2),
        "Rosetta_total_score": rosetta_output['total_score'],
        "Rosetta_interaction_score": rosetta_output['interaction_score'],
        "Rosetta_interaction_score_precomputed": rosetta_output['rosetta_computed'],
    }
    # write to interface info file
    save_json(interfaceInfo, interface_info_file)
    with open(interface_residue_file, "w") as out:
        header = "ChainID,IntID,ResID,Type,Potential,SASA,RelSASA,IsHotspot,InsCode"
        out.write(header+"\n")
        for id in interface_residues1:
            residue = residues1.get_residue_by_id(id)
            c = residue.chain_id
            iid = residue.id
            resid = residue.res_id
            aa = residue.res_name
            pot = residue.pair_potential
            sasa = "NA" if residue.sasa is None else residue.sasa
            relsasa = "NA" if residue.relative_sasa is None else residue.relative_sasa
            if None in [c, iid, resid, aa, pot]:
                continue
            is_hotspot = "True" if residue.is_hotspot else "False"
            inscode = "None" if residue.ins_code is None else residue.ins_code
            out.write("{},{},{},{},{},{},{},{},{}\n".format(
                c,iid,resid,aa,pot,sasa,relsasa,is_hotspot,inscode
            ))
        for id in interface_residues2:
            residue = residues2.get_residue_by_id(id)
            c = residue.chain_id
            iid = residue.id
            resid = residue.res_id
            aa = residue.res_name
            pot = residue.pair_potential
            sasa = "NA" if residue.sasa is None else residue.sasa
            relsasa = "NA" if residue.relative_sasa is None else residue.relative_sasa
            if None in [c, iid, resid, aa, pot]:
                continue
            is_hotspot = "True" if residue.is_hotspot else "False"
            inscode = "None" if residue.ins_code is None else residue.ins_code
            out.write("{},{},{},{},{},{},{},{},{}\n".format(
                c,iid,resid,aa,pot,sasa,relsasa,is_hotspot,inscode
            ))
    filenames = {
        'twochain_file': twochain_file,
        'chain1_file': chain1_file,
        'chain2_file': chain2_file,
        'intres_file': interface_residue_file,
        'info_file': interface_info_file,
    }
    return filenames

def parse_args():
    """Parse arguments."""
    parser = argparse.ArgumentParser(
      description="Prepare template/interface/hotspot for single PDB entry."
    )
    parser.add_argument(
      '--pdb', type=str, required=False,
      help="A single 4-letter PDB id to process."
    )
    parser.add_argument(
        '--mode', type=str, default="all",
        help="Default --mode=all, parse all assemblies. Set to dataprep and provide pdb file. Set to anything else to parse only the given chains."
    )
    parser.add_argument(
        '--pdbfile', type=str, required=False,
        default='',
        help="Provide .pdb file when --mode=dataprep. Intended to prepare interfaces from existing databases (e.g. dockground)"
    )
    parser.add_argument(
        '--chain1', type=str,
        help="chain id 1. Must be specified if --mode is not all"
    )
    parser.add_argument(
        '--chain2', type=str,
        help='chain id 2. Must be specified if --mode is not all'
    )

    parser.add_argument(
      '--pdbpath', type=str,
      default="/data/limh5/pdb",
      help="Path to download CIF."
    )
    parser.add_argument(
      '--templatepath', type=str,
      default="/data/limh5/hmi_templates_v4/interfaces",
      help="Base path to processed interfaces."
    )
    parser.add_argument(
        '--log', default="INFO", help="Logging level. Set to DEBUG for more details."
    )
    return parser.parse_args()


def main(args):
    FORMAT = '%(asctime)-15s %(message)s'
    logging.basicConfig(format=FORMAT, level=getattr(logging, args.log.upper()))

    # --- directly from settings
    ROOT_DIR = Path(__file__).resolve(strict=True).parent.parent.parent
    APPS_DIR = ROOT_DIR / "hmi_pred"
    MEDIA_ROOT = str(APPS_DIR / "media")
    
    #ROSETTA_BASEPATH = Path("/data/limh5/rosetta_bin_linux_2021.16.61629_bundle").resolve()
    #ROSETTA_BASEPATH = Path("/Users/limh5/rosetta_bin_mac_2021.16.61629_bundle").resolve()
    ROSETTA_BASEPATH = Path("/home/hansaim/rosetta/rosetta_bin_linux_2021.16.61629_bundle").resolve()
    EPPIC_BASEPATH = Path("/data/limh5/eppic-cli-3.2.5/").resolve()
    HMI_TOOLS = {
    # PATHS TO EXECUTABLE TOOLS AND RELEVANT FILES
        'EPPIC': os.path.join(EPPIC_BASEPATH, 'bin/eppic'),
        'EPPIC_CONF': os.path.join(EPPIC_BASEPATH, '.eppic.conf'),
        'TMALIGN': os.path.join(ROOT_DIR, 'scripts/bin/TMalign'),
        'ROSETTAPREPACK': os.path.join(ROSETTA_BASEPATH, 'main/source/bin/docking_prepack_protocol.static.linuxgccrelease'),
        'ROSETTADOCK': os.path.join(ROSETTA_BASEPATH, 'main/source/bin/docking_protocol.static.linuxgccrelease'),
        #'ROSETTAPREPACK': os.path.join(ROSETTA_BASEPATH, 'main/source/bin/docking_prepack_protocol.static.macosclangrelease'),
        #'ROSETTADOCK': os.path.join(ROSETTA_BASEPATH, 'main/source/bin/docking_protocol.static.macosclangrelease'),
        'ROSETTA_DB': os.path.join(ROSETTA_BASEPATH, 'main/database/'),
    }
    # --- directly from settings

    if args.mode.lower() in ['dataprep', 'data prep', 'database']:
        # collect single interface from .pdb file in a database
        # args.pdbfile must be given
        # no need to parse gene info
        logging.info("Dataprep mode.")
        interfaceDir = os.path.dirname(args.pdbfile)
        pdbid = args.pdbfile.strip().split('/')[-1].split('.')[0].replace('_st','')
        chain1 = args.chain1
        chain2 = args.chain2
        exists, ne = check_interface_files(interfaceDir, pdbid, chain1, chain2)
        if exists:
            # all files for the given interface exists; skip
            print("Interface exists for {}_{}_{} in directory: {}".format(pdbid, chain1, chain2, interfaceDir))
        else:
            logging.info("Collecting interface from {}".format(args.pdbfile))
            is_valid, interface = collect_interfaces_from_cif(
                                        args.pdbfile, assembly_id=None, chain_mode='dataprep', chain1=chain1, chain2=chain2,
                                        sasa_probe_radius=1.4, sasa_point_number=1000, sasa_ignore_ions=True,
                                        scaffold_distance_threshold=6.0,
                                        hotspot_com_distance_threshold=7.0,
                                        hotspot_res_id_diff=4,
                                        hotspot_min_total_potential=18.0,
                                        hotspot_max_sasa=20.0,
                                        hotregion_com_distance_threshold=6.5,
                                        combined_file_output_dir=interfaceDir,
                                        clean=False,)

            _info = ('','','','','') # dummy info
            chain_info = {chain1:_info, chain2:_info}

            if is_valid:
                logging.info("Interface is valid.")
            else:
                logging.info("Interface invalid.")
            logging.info("Writing interface files to {}".format(interfaceDir))
            filenames = write_interface_files(interface, interfaceDir, pdbid, (chain1, chain2), chain_info, 
                assembly_id=None, HMI_TOOLS=HMI_TOOLS, validate=False) 
            logging.info("Interface extraction complete. Written to {}".format(interfaceDir))
            logging.info("    Filenames: {}".format(filenames))
    else:

        pdbid = args.pdb.upper()
        pdbBasePath = Path(args.pdbpath).resolve()
        pdbDir = os.path.join(pdbBasePath, pdbid[1:3])
        if not Path(pdbDir).is_dir():
            Popen("mkdir -p {}".format(pdbDir), shell=True).wait()

        cif_file = os.path.join(pdbDir, "{}.cif".format(pdbid))
        
        if not check_downloaded_cif(cif_file):
            # cif file does not exist: not pre-downloaded
            logging.debug("{} does not exist. Try fresh download.".format(cif_file))
            cif_file = utils.download_pdb(pdbid, pdbDir, cif=True, decompress=True)
            if cif_file is None:
                # not pre-downloaded, fresh download unsuccessful
                raise ValueError("Download unsuccessful.")
        else:
            logging.debug("{} found. Continuing with this file.".format(cif_file))
            
        templateBasePath = Path(args.templatepath).resolve()
        interfaceDir = os.path.join(templateBasePath, pdbid[1:3])
        if not Path(interfaceDir).is_dir():
            Popen("mkdir -p {}".format(interfaceDir), shell=True).wait()

        chain_info = protein.get_chain_info(cif_file)
        chain_ids = sorted(list(chain_info.keys()))
        # {'A': ('P20231', 'TRYB2_HUMAN', '9606', 'Homo sapiens', 'human'),
        # 'B': ('P20231', 'TRYB2_HUMAN', '9606', 'Homo sapiens', 'human'),
        # 'C': ('P20231', 'TRYB2_HUMAN', '9606', 'Homo sapiens', 'human'),
        # 'D': ('P20231', 'TRYB2_HUMAN', '9606', 'Homo sapiens', 'human')}

        logging.debug("PDB: {}. Chain info: {}".format(pdbid, chain_info))

        chain_pairs = []  # at least one human chain
        for i in range(len(chain_ids)):
            c1 = chain_ids[i]
            c1_taxid = chain_info[c1][2]

            for j in range(len(chain_ids)):
                if i >= j:
                    continue
                c2 = chain_ids[j]
                c2_taxid = chain_info[c2][2]
                tax_ids = [c1_taxid, c2_taxid]
                if '9606' not in tax_ids:
                    # non-human chain combination; skip
                    continue
                else:
                    # at least one human chain
                    chain_pairs.append((c1,c2))

        obj = pdbx.PDBxFile.read(cif_file)
        try:
            assembly_list = pdbx.list_assemblies(obj)
            assemblies = {k:pdbx.get_assembly(obj, k, altloc='occupancy') for k in assembly_list}
            assem_chains = {k:list(set(assemblies[k][0].chain_id)) for k in assemblies}

            chainpair2assemblyid = {}
            for pair in chain_pairs:
                for aid in assem_chains:
                    if (pair[0] in assem_chains[aid]) and (pair[1] in assem_chains[aid]):
                        chainpair2assemblyid[pair] = aid
                        break
            assembly_available = True
        except:
            assembly_available = False

        for pair in chain_pairs:
            if args.mode.lower() not in ['all', 'all_assemblies']:
                pair_ = tuple(sorted([args.chain1, args.chain2]))
                if (pair[0]!=pair_[0]) or (pair[1]!=pair_[1]):
                    continue
            if not assembly_available:
                assembly_id = None
            else:
                if pair not in chainpair2assemblyid:
                    logging.debug("Chain pair {} could not be found in any assembly in {}".format(pair, cif_file))
                    continue
                assembly_id = chainpair2assemblyid[pair]
            is_valid, interface = collect_interfaces_from_cif(
                                    cif_file, assembly_id=assembly_id, chain_mode='specific', chain1=pair[0], chain2=pair[1],
                                    sasa_probe_radius=1.4, sasa_point_number=1000, sasa_ignore_ions=True,
                                    scaffold_distance_threshold=6.0,
                                    hotspot_com_distance_threshold=7.0,
                                    hotspot_res_id_diff=4,
                                    hotspot_min_total_potential=18.0,
                                    hotspot_max_sasa=20.0,
                                    hotregion_com_distance_threshold=6.5,
                                    combined_file_output_dir=interfaceDir,
                                    clean=False,
                                    
                                )
            if not is_valid:
                print(interface) #error message
                continue
            filenames = write_interface_files(interface, interfaceDir, pdbid, pair, chain_info, assembly_id=assembly_id, HMI_TOOLS=HMI_TOOLS) 
    

if __name__ == "__main__":
    args = parse_args()
    main(args)
