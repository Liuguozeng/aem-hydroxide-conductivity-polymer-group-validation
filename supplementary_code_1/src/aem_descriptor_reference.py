"""Reference implementation of the 78 computed AEM descriptors.

The function accepts one or two repeat-unit SMILES strings and their composition
weights. No literature-derived structures or observation records are included.
"""
from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

PATTERNS = {'ether': Chem.MolFromSmarts('[#6]-O-[#6]'), 'aryl_ether': Chem.MolFromSmarts('a-O-a'), 'sulfone': Chem.MolFromSmarts('[#16X4](=[OX1])(=[OX1])'), 'ketone': Chem.MolFromSmarts('[#6][CX3](=O)[#6]'), 'carbonyl': Chem.MolFromSmarts('[CX3]=[OX1]'), 'hydroxyl': Chem.MolFromSmarts('[OX2H]'), 'guanidinium': Chem.MolFromSmarts('[NX3][CX3](=[NX3+])[NX3]'), 'benzylic_link': Chem.MolFromSmarts('[c]-[CH2,CH]-[N+,P+,S+]')}

def to_float_or_nan(value: Any) -> float:
    if value is None:
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    text = str(value).strip().replace('%', '')
    if text == '' or text.lower() in {'nan', 'none'}:
        return np.nan
    try:
        return float(text)
    except Exception:
        return np.nan

def normalize_ratios(ratio_1: Any, ratio_2: Any, repeating_unit_1: str, repeating_unit_2: str) -> tuple[float, float, str]:
    ru1 = '' if repeating_unit_1 is None else str(repeating_unit_1).strip()
    ru2 = '' if repeating_unit_2 is None else str(repeating_unit_2).strip()
    r1 = to_float_or_nan(ratio_1)
    r2 = to_float_or_nan(ratio_2)
    if ru1 == ru2:
        return (1.0, 0.0, 'homopolymer_same_unit')
    if np.isnan(r1) and np.isnan(r2):
        if ru2 != '':
            return (0.5, 0.5, 'both_missing_imputed_50_50')
        return (1.0, 0.0, 'single_unit_default')
    if np.isnan(r1) and (not np.isnan(r2)):
        r1 = max(0.0, 1.0 - r2)
        return (float(r1), float(r2), 'ratio_1_filled_from_ratio_2')
    if np.isnan(r2) and (not np.isnan(r1)):
        r2 = max(0.0, 1.0 - r1)
        return (float(r1), float(r2), 'ratio_2_filled_from_ratio_1')
    total = r1 + r2
    if total <= 0:
        return (0.5, 0.5, 'invalid_total_reset_50_50')
    return (float(r1 / total), float(r2 / total), 'given')

def canonical_smiles(smiles: Any) -> str:
    if smiles is None or (isinstance(smiles, float) and math.isnan(smiles)):
        return ''
    text = str(smiles).strip()
    if text == '' or text.lower() == 'nan':
        return ''
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return text
    return Chem.MolToSmiles(mol)

def _mol_from_smiles(smiles: str) -> Chem.Mol | None:
    if smiles is None:
        return None
    text = str(smiles).strip()
    if text == '' or text.lower() == 'nan':
        return None
    mol = Chem.MolFromSmiles(text)
    return mol

def _count_matches(mol: Chem.Mol, pattern_name: str) -> int:
    pattern = PATTERNS[pattern_name]
    return len(mol.GetSubstructMatches(pattern)) if pattern is not None else 0

def _atom_fraction(mol: Chem.Mol, predicate) -> float:
    atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 0]
    if not atoms:
        return np.nan
    return sum((1 for atom in atoms if predicate(atom))) / len(atoms)

def _shortest_distance_to_targets(mol: Chem.Mol, start_idx: int, target_indices: Sequence[int]) -> float:
    if not target_indices:
        return np.nan
    distances = []
    for target_idx in target_indices:
        if target_idx == start_idx:
            continue
        try:
            path = Chem.rdmolops.GetShortestPath(mol, start_idx, target_idx)
        except Exception:
            continue
        if path:
            distances.append(len(path) - 1)
    return float(min(distances)) if distances else np.nan

def _weighted_scalar_merge(a: dict[str, float], b: dict[str, float], w_a: float, w_b: float) -> dict[str, float]:
    keys = sorted(set(a) | set(b))
    merged = {}
    for key in keys:
        va = a.get(key, np.nan)
        vb = b.get(key, np.nan)
        if pd.isna(va) and pd.isna(vb):
            merged[key] = np.nan
        elif pd.isna(va):
            merged[key] = vb
        elif pd.isna(vb):
            merged[key] = va
        else:
            merged[key] = w_a * va + w_b * vb
    return merged

def _is_double_bonded_hetero(atom: Chem.Atom) -> bool:
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if bond.GetBondTypeAsDouble() > 1.1 and other.GetAtomicNum() in {6, 7, 8, 16}:
            return True
    return False

def _cation_family_flags(mol: Chem.Mol) -> dict[str, float]:
    ring_info = mol.GetRingInfo().AtomRings()
    flags = {'mech__has_quaternary_ammonium_like': 0.0, 'mech__has_imidazolium_like': 0.0, 'mech__has_piperidinium_like': 0.0, 'mech__has_pyrrolidinium_like': 0.0, 'mech__has_pyridinium_like': 0.0, 'mech__has_phosphonium_like': 0.0, 'mech__has_guanidinium_like': float(_count_matches(mol, 'guanidinium') > 0)}
    for atom in _candidate_cation_atoms(mol):
        atomic_num = atom.GetAtomicNum()
        atom_rings = [ring for ring in ring_info if atom.GetIdx() in ring]
        ring_sizes = [len(ring) for ring in atom_rings]
        ring_has_second_n = False
        for ring in atom_rings:
            n_count = sum((1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7))
            if n_count >= 2:
                ring_has_second_n = True
                break
        if atomic_num == 15:
            flags['mech__has_phosphonium_like'] = 1.0
        elif atomic_num == 7:
            if not ring_sizes:
                flags['mech__has_quaternary_ammonium_like'] = 1.0
            else:
                min_ring = min(ring_sizes)
                if atom.GetIsAromatic() and min_ring == 6:
                    flags['mech__has_pyridinium_like'] = 1.0
                elif ring_has_second_n and min_ring == 5:
                    flags['mech__has_imidazolium_like'] = 1.0
                elif min_ring == 6 and (not atom.GetIsAromatic()):
                    flags['mech__has_piperidinium_like'] = 1.0
                elif min_ring == 5 and (not atom.GetIsAromatic()):
                    flags['mech__has_pyrrolidinium_like'] = 1.0
                else:
                    flags['mech__has_quaternary_ammonium_like'] = 1.0
    flags['mech__cation_family_count'] = float(sum((v > 0 for v in flags.values())))
    return flags

def _pairwise_shortest_distances(mol: Chem.Mol, atom_indices: Sequence[int]) -> list[float]:
    dists = []
    atom_indices = list(dict.fromkeys((int(i) for i in atom_indices)))
    for i in range(len(atom_indices)):
        for j in range(i + 1, len(atom_indices)):
            try:
                path = Chem.rdmolops.GetShortestPath(mol, atom_indices[i], atom_indices[j])
            except Exception:
                path = ()
            if path:
                dists.append(float(len(path) - 1))
    return dists

def _cation_environment_topology(mol: Chem.Mol, cation_atoms: Sequence[Chem.Atom]) -> dict[str, float]:
    if not cation_atoms:
        return {'mech__cation_is_pendant_fraction': 0.0, 'mech__cation_is_backbone_bound_fraction': 0.0, 'mech__cation_to_anchor_dist_mean': 0.0, 'mech__cation_to_anchor_dist_min': 0.0, 'mech__cation_to_aromatic_dist_mean': 0.0, 'mech__cation_to_cation_dist_mean': 0.0, 'mech__spacer_oeg_unit_count': 0.0, 'mech__spacer_branching_index': 0.0, 'mech__beta_h_accessibility_proxy': 0.0, 'mech__ring_cation_fraction': 0.0}
    dummy_targets = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    aromatic_targets = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    anchor_dists = []
    aromatic_dists = []
    pendant_flags = []
    backbone_flags = []
    ring_flags = []
    oeg_counts = []
    branching_scores = []
    beta_h_scores = []
    for atom in cation_atoms:
        idx = atom.GetIdx()
        ring_flags.append(float(atom.IsInRing()))
        anchor_targets = dummy_targets if dummy_targets else aromatic_targets
        anchor_dist = _shortest_distance_to_targets(mol, idx, anchor_targets)
        aromatic_dist = _shortest_distance_to_targets(mol, idx, aromatic_targets)
        if np.isnan(anchor_dist):
            anchor_dist = aromatic_dist if not np.isnan(aromatic_dist) else 0.0
        if np.isnan(aromatic_dist):
            aromatic_dist = anchor_dist if not np.isnan(anchor_dist) else 0.0
        anchor_dists.append(float(anchor_dist))
        aromatic_dists.append(float(aromatic_dist))
        pendant_flags.append(float(anchor_dist >= 2.0))
        backbone_flags.append(float(anchor_dist <= 1.0 or aromatic_dist <= 1.0))
        visited = {idx}
        frontier = {idx}
        hetero_oxygen_count = 0
        branch_count = 0
        beta_h = 0.0
        for depth in range(1, 4):
            next_frontier = set()
            for node_idx in frontier:
                node = mol.GetAtomWithIdx(node_idx)
                for nb in node.GetNeighbors():
                    nb_idx = nb.GetIdx()
                    if nb_idx in visited:
                        continue
                    visited.add(nb_idx)
                    next_frontier.add(nb_idx)
                    if nb.GetAtomicNum() == 8:
                        hetero_oxygen_count += 1
                    if nb.GetAtomicNum() == 6 and nb.GetDegree() >= 3:
                        branch_count += 1
                    if depth <= 2 and nb.GetAtomicNum() == 6:
                        if not nb.GetIsAromatic() and nb.GetTotalNumHs() > 0:
                            beta_h += 1.0
            frontier = next_frontier
        oeg_counts.append(float(hetero_oxygen_count))
        branching_scores.append(float(branch_count))
        beta_h_scores.append(float(beta_h > 0))
    cation_distances = _pairwise_shortest_distances(mol, [atom.GetIdx() for atom in cation_atoms])
    return {'mech__cation_is_pendant_fraction': float(np.mean(pendant_flags)), 'mech__cation_is_backbone_bound_fraction': float(np.mean(backbone_flags)), 'mech__cation_to_anchor_dist_mean': float(np.mean(anchor_dists)), 'mech__cation_to_anchor_dist_min': float(np.min(anchor_dists)), 'mech__cation_to_aromatic_dist_mean': float(np.mean(aromatic_dists)), 'mech__cation_to_cation_dist_mean': float(np.mean(cation_distances)) if cation_distances else 0.0, 'mech__spacer_oeg_unit_count': float(np.mean(oeg_counts)), 'mech__spacer_branching_index': float(np.mean(branching_scores)), 'mech__beta_h_accessibility_proxy': float(np.mean(beta_h_scores)), 'mech__ring_cation_fraction': float(np.mean(ring_flags))}

def _mechanism_core_descriptors(mol: Chem.Mol) -> dict[str, float]:
    heavy_atoms = max(float(mol.GetNumHeavyAtoms()), 1.0)
    hetero_count = float(sum((1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (0, 1, 6))))
    cation_atoms = _candidate_cation_atoms(mol)
    ionogenic_atoms = _candidate_ionogenic_site_atoms(mol)
    site_atoms = ionogenic_atoms if ionogenic_atoms else cation_atoms
    n_cation = float(len(cation_atoms))
    n_ionogenic = float(len(ionogenic_atoms))
    effective_site_count = max(n_cation, n_ionogenic, 1.0)
    n_site_safe = max(effective_site_count, 1.0)
    ether_count = float(_count_matches(mol, 'ether'))
    aryl_ether_count = float(_count_matches(mol, 'aryl_ether'))
    sulfone_count = float(_count_matches(mol, 'sulfone'))
    ketone_count = float(_count_matches(mol, 'ketone'))
    carbonyl_count = float(_count_matches(mol, 'carbonyl'))
    benzylic_proxy = float(_count_matches(mol, 'benzylic_link'))
    hydroxyl_count = float(_count_matches(mol, 'hydroxyl'))
    aromatic_fraction = _atom_fraction(mol, lambda atom: atom.GetIsAromatic())
    ring_count = float(rdMolDescriptors.CalcNumRings(mol))
    aromatic_rings = float(rdMolDescriptors.CalcNumAromaticRings(mol))
    rotatable_bonds = float(Lipinski.NumRotatableBonds(mol))
    logp = float(Descriptors.MolLogP(mol))
    fluorine_fraction = _atom_fraction(mol, lambda atom: atom.GetAtomicNum() == 9)
    frac_csp3 = float(rdMolDescriptors.CalcFractionCSP3(mol))
    spiro_atoms = float(getattr(rdMolDescriptors, 'CalcNumSpiroAtoms')(mol))
    bridgehead_atoms = float(getattr(rdMolDescriptors, 'CalcNumBridgeheadAtoms')(mol))
    topo = _cation_environment_topology(mol, site_atoms)
    cation_density = effective_site_count / heavy_atoms
    hydrophilic_density = (hetero_count + 1.6 * effective_site_count + hydroxyl_count + topo['mech__spacer_oeg_unit_count']) / heavy_atoms
    hydrophobicity_proxy = (max(logp, 0.0) + aromatic_fraction + max(fluorine_fraction, 0.0) + (1.0 - min(frac_csp3, 1.0))) / 4.0
    flexibility_proxy = (rotatable_bonds + topo['mech__spacer_oeg_unit_count']) / heavy_atoms
    rigidity_proxy = (aromatic_rings + ring_count + ketone_count + sulfone_count + 1.0) / (rotatable_bonds + 1.0)
    steric_hindrance_proxy = topo['mech__ring_cation_fraction'] + topo['mech__spacer_branching_index'] / max(heavy_atoms, 1.0)
    sn2_risk_proxy = benzylic_proxy / n_site_safe + aryl_ether_count / heavy_atoms + max(0.0, 2.0 - topo['mech__cation_to_anchor_dist_mean']) / 4.0
    hofmann_e2_risk_proxy = topo['mech__beta_h_accessibility_proxy'] * (1.0 - 0.5 * topo['mech__ring_cation_fraction']) + 0.25 * (1.0 - topo['mech__cation_is_pendant_fraction'])
    free_volume_proxy = (spiro_atoms + bridgehead_atoms + fluorine_fraction * heavy_atoms + 0.5 * frac_csp3 * heavy_atoms) / (heavy_atoms + 1.0)
    alkaline_vulnerability_proxy = 0.45 * sn2_risk_proxy + 0.35 * hofmann_e2_risk_proxy + 0.2 * carbonyl_count / heavy_atoms + 0.25 * float(_count_matches(mol, 'guanidinium') > 0)
    pendant_fraction = topo['mech__cation_is_pendant_fraction']
    backbone_bound_fraction = topo['mech__cation_is_backbone_bound_fraction']
    ring_cation_fraction = topo['mech__ring_cation_fraction']
    spacer_oeg = topo['mech__spacer_oeg_unit_count']
    spacer_len = topo['mech__cation_to_anchor_dist_mean']
    neutral_ionogenic_fraction = max(0.0, n_ionogenic - n_cation) / max(n_ionogenic, 1.0)
    cation_accessibility_proxy = (0.4 * pendant_fraction + 0.2 * min(spacer_len / 6.0, 1.5) + 0.18 * min(flexibility_proxy * 3.0, 1.5) + 0.12 * min((spacer_oeg + 1.0) / 3.0, 1.5) + 0.1 * (1.0 - neutral_ionogenic_fraction)) / (1.0 + 0.6 * steric_hindrance_proxy + 0.45 * ring_cation_fraction + 0.35 * backbone_bound_fraction)
    water_network_proxy = hydrophilic_density * (1.0 + 0.35 * spacer_oeg + 0.2 * pendant_fraction + 0.1 * neutral_ionogenic_fraction) / (1.0 + 0.55 * hydrophobicity_proxy + 0.2 * backbone_bound_fraction)
    cation_localization_penalty = backbone_bound_fraction + 0.5 * ring_cation_fraction + 0.25 * steric_hindrance_proxy
    hydration_shell_overlap_proxy = water_network_proxy * (1.0 + 0.45 * pendant_fraction + 0.3 * min(spacer_len / 6.0, 1.5) + 0.2 * min((spacer_oeg + 1.0) / 2.5, 1.5)) / (1.0 + 0.3 * hydrophobicity_proxy + 0.35 * cation_localization_penalty)
    grotthuss_window_proxy = hydration_shell_overlap_proxy * (1.0 + 0.3 * hydrophilic_density + 0.2 * flexibility_proxy) / (1.0 + 0.35 * steric_hindrance_proxy + 0.25 * hydrophobicity_proxy)
    channel_index = (0.28 * free_volume_proxy + 0.24 * water_network_proxy + 0.18 * cation_accessibility_proxy + 0.16 * hydration_shell_overlap_proxy + 0.14 * grotthuss_window_proxy) * (1.0 + 0.18 * flexibility_proxy) / (1.0 + 0.35 * cation_localization_penalty)
    transport_stability_balance = (0.65 * channel_index + 0.35 * grotthuss_window_proxy) / (1.0 + alkaline_vulnerability_proxy)
    mech = {'mech__cation_count': n_cation, 'mech__ionogenic_site_count_proxy': n_ionogenic, 'mech__effective_ion_site_density': cation_density, 'mech__neutral_ionogenic_fraction': neutral_ionogenic_fraction, 'mech__multication_repeat_proxy': float(effective_site_count >= 2), 'mech__local_charge_density_proxy': cation_density, 'mech__cyclic_cation_fraction': ring_cation_fraction, 'mech__benzylic_link_density': benzylic_proxy / n_site_safe, 'mech__spacer_length_proxy': spacer_len, 'mech__cation_hetero_neighborhood_proxy': spacer_oeg / max(heavy_atoms, 1.0), 'mech__steric_hindrance_proxy': steric_hindrance_proxy, 'mech__aromaticity_proxy': aromatic_fraction, 'mech__backbone_rigidity_proxy': rigidity_proxy, 'mech__sidechain_flexibility_proxy': flexibility_proxy, 'mech__ether_motif_density': ether_count / heavy_atoms, 'mech__aryl_ether_motif_density': aryl_ether_count / heavy_atoms, 'mech__ether_free_backbone_proxy': float(aryl_ether_count == 0 and ether_count == 0), 'mech__sulfone_motif_density': sulfone_count / heavy_atoms, 'mech__ketone_motif_density': ketone_count / heavy_atoms, 'mech__hydrophobicity_proxy': hydrophobicity_proxy, 'mech__fluorination_proxy': fluorine_fraction, 'mech__hydrophilic_density_proxy': hydrophilic_density, 'mech__cation_accessibility_proxy': cation_accessibility_proxy, 'mech__water_network_proxy': water_network_proxy, 'mech__hydration_shell_overlap_proxy': hydration_shell_overlap_proxy, 'mech__grotthuss_window_proxy': grotthuss_window_proxy, 'mech__channel_index': channel_index, 'mech__transport_stability_balance': transport_stability_balance, 'mech__cation_localization_penalty': cation_localization_penalty, 'mech__sn2_risk_proxy': sn2_risk_proxy, 'mech__hofmann_e2_risk_proxy': hofmann_e2_risk_proxy, 'mech__free_volume_proxy': free_volume_proxy, 'mech__alkaline_vulnerability_proxy': alkaline_vulnerability_proxy}
    mech.update(topo)
    mech.update(_cation_family_flags(mol))
    return mech

def _weighted_harmonic_value(value_1: float, value_2: float, weight_1: float, weight_2: float) -> float:
    if pd.isna(value_1) and pd.isna(value_2):
        return np.nan
    pairs: list[tuple[float, float]] = []
    if pd.notna(value_1) and float(value_1) > 0:
        pairs.append((float(value_1), float(weight_1)))
    if pd.notna(value_2) and float(value_2) > 0:
        pairs.append((float(value_2), float(weight_2)))
    if not pairs:
        return np.nan
    if len(pairs) == 1:
        return pairs[0][0]
    denom = sum((w / max(v, 1e-12) for v, w in pairs))
    if denom <= 0:
        return np.nan
    return float(sum((w for _, w in pairs)) / denom)

def _composition_aware_merge(desc1: dict[str, float], desc2: dict[str, float], r1: float, r2: float) -> dict[str, float]:
    scalar_keys = [k for k in desc1 if not (k.startswith('fpb__') or k.startswith('fpc__'))]
    fpb_keys = [k for k in desc1 if k.startswith('fpb__')]
    fpc_keys = [k for k in desc1 if k.startswith('fpc__')]
    out = _weighted_scalar_merge({k: desc1[k] for k in scalar_keys}, {k: desc2[k] for k in scalar_keys}, r1, r2)
    contrast_pairs = {'mix__cation_localization_contrast': 'mech__local_charge_density_proxy', 'mix__hydrophilic_contrast': 'mech__hydrophilic_density_proxy', 'mix__hydrophobicity_contrast': 'mech__hydrophobicity_proxy', 'mix__rigidity_contrast': 'mech__backbone_rigidity_proxy', 'mix__spacer_contrast': 'mech__spacer_length_proxy', 'mix__fluorination_contrast': 'mech__fluorination_proxy', 'mix__free_volume_contrast': 'mech__free_volume_proxy', 'mix__alkaline_risk_contrast': 'mech__alkaline_vulnerability_proxy', 'mix__accessibility_contrast': 'mech__cation_accessibility_proxy', 'mix__water_network_contrast': 'mech__water_network_proxy', 'mix__hydration_overlap_contrast': 'mech__hydration_shell_overlap_proxy', 'mix__grotthuss_window_contrast': 'mech__grotthuss_window_proxy', 'mix__channel_index_contrast': 'mech__channel_index', 'mix__transport_balance_contrast': 'mech__transport_stability_balance', 'mix__localization_penalty_contrast': 'mech__cation_localization_penalty'}
    for out_name, key in contrast_pairs.items():
        v1 = desc1.get(key, np.nan)
        v2 = desc2.get(key, np.nan)
        out[out_name] = np.nan if pd.isna(v1) or pd.isna(v2) else abs(float(v1) - float(v2))
    for key in fpb_keys:
        idx = key.split('__', 1)[1]
        v1 = 0.0 if pd.isna(desc1[key]) else float(desc1[key])
        v2 = 0.0 if pd.isna(desc2[key]) else float(desc2[key])
        out[f'fpm__{idx}'] = r1 * v1 + r2 * v2
        out[f'fpx__{idx}'] = abs(v1 - v2)
    for key in fpc_keys:
        idx = key.split('__', 1)[1]
        v1 = 0.0 if pd.isna(desc1[key]) else float(desc1[key])
        v2 = 0.0 if pd.isna(desc2[key]) else float(desc2[key])
        out[f'fpcm__{idx}'] = r1 * v1 + r2 * v2
    is_copolymer = float(min(r1, r2) > 0)
    composition_balance = 1.0 - abs(r1 - r2)
    out['mix__asymmetric_copolymer_flag'] = float(is_copolymer > 0 and abs(r1 - r2) >= 0.3)
    out['mix__ionic_site_localization_flag'] = float(out.get('mix__cation_localization_contrast', 0.0) >= 0.08 and is_copolymer > 0)
    water_h = _weighted_harmonic_value(desc1.get('mech__water_network_proxy', np.nan), desc2.get('mech__water_network_proxy', np.nan), r1, r2)
    channel_h = _weighted_harmonic_value(desc1.get('mech__channel_index', np.nan), desc2.get('mech__channel_index', np.nan), r1, r2)
    overlap_h = _weighted_harmonic_value(desc1.get('mech__hydration_shell_overlap_proxy', np.nan), desc2.get('mech__hydration_shell_overlap_proxy', np.nan), r1, r2)
    out['mix__water_network_harmonic'] = water_h
    out['mix__channel_index_harmonic'] = channel_h
    out['mix__hydration_overlap_harmonic'] = overlap_h
    weighted_channel = r1 * float(desc1.get('mech__channel_index', 0.0) or 0.0) + r2 * float(desc2.get('mech__channel_index', 0.0) or 0.0)
    out['mix__channel_bottleneck_penalty'] = max(0.0, weighted_channel - (0.0 if pd.isna(channel_h) else float(channel_h)))
    contrast_terms = [out.get('mix__hydrophilic_contrast', 0.0), out.get('mix__hydrophobicity_contrast', 0.0), out.get('mix__rigidity_contrast', 0.0), out.get('mix__free_volume_contrast', 0.0)]
    out['mix__segregation_index'] = composition_balance * float(np.nanmean(contrast_terms))
    heterogeneity_terms = [out.get('mix__hydrophilic_contrast', 0.0), out.get('mix__accessibility_contrast', 0.0), out.get('mix__channel_index_contrast', 0.0), out.get('mix__localization_penalty_contrast', 0.0), out.get('mix__channel_bottleneck_penalty', 0.0)]
    out['mix__transport_heterogeneity'] = float(np.nanmean(heterogeneity_terms))
    channel_terms = [out.get('mix__channel_index_harmonic', np.nan), out.get('mix__water_network_harmonic', np.nan), out.get('mix__hydration_overlap_harmonic', np.nan), out.get('mix__free_volume_contrast', 0.0), out.get('mix__accessibility_contrast', 0.0)]
    out['mix__channel_formation_index'] = composition_balance * float(np.nanmean(channel_terms))
    hydration_terms = [out.get('mix__water_network_harmonic', np.nan), out.get('mix__hydration_overlap_harmonic', np.nan), out.get('mix__grotthuss_window_contrast', 0.0)]
    out['mix__hydration_overlap_index'] = composition_balance * float(np.nanmean(hydration_terms))
    arch_weight = 1.15 if is_copolymer > 0 else 1.0
    out['mix__architecture_weighted_segregation'] = arch_weight * out['mix__segregation_index']
    out['mix__architecture_weighted_channel_index'] = arch_weight * out['mix__channel_formation_index']
    return out

def _candidate_cation_atoms(mol: Chem.Mol) -> list[Chem.Atom]:
    ring_info = mol.GetRingInfo().AtomRings()
    atoms: list[Chem.Atom] = []
    seen: set[int] = set()

    def add(atom: Chem.Atom):
        idx = atom.GetIdx()
        if idx not in seen:
            seen.add(idx)
            atoms.append(atom)
    for atom in mol.GetAtoms():
        anum = atom.GetAtomicNum()
        if atom.GetFormalCharge() > 0:
            add(atom)
            continue
        if anum == 15 and atom.GetTotalDegree() >= 4:
            add(atom)
            continue
        if anum != 7:
            continue
        neighbors = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1]
        carbon_neighbors = sum((nb.GetAtomicNum() == 6 for nb in neighbors))
        exocyclic_carbons = sum((nb.GetAtomicNum() == 6 and (not nb.IsInRing()) for nb in neighbors))
        atom_rings = [ring for ring in ring_info if atom.GetIdx() in ring]
        ring_n_counts = [sum((mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in ring)) for ring in atom_rings]
        aromatic_5n = any((len(ring) == 5 and n_count >= 2 for ring, n_count in zip(atom_rings, ring_n_counts)))
        aza_fused = any((len(ring) in (5, 6) and n_count >= 2 for ring, n_count in zip(atom_rings, ring_n_counts)))
        sat6n = any((len(ring) == 6 and (not any((mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))) for ring in atom_rings))
        sat5n = any((len(ring) == 5 and (not any((mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))) for ring in atom_rings))
        if aromatic_5n or aza_fused:
            add(atom)
            continue
        if _is_double_bonded_hetero(atom):
            continue
        if atom.GetTotalDegree() >= 4:
            add(atom)
            continue
        if atom.GetTotalDegree() == 3 and carbon_neighbors >= 2:
            add(atom)
            continue
        if sat6n or sat5n:
            add(atom)
            continue
        if aza_fused and (carbon_neighbors >= 1 or exocyclic_carbons >= 1):
            add(atom)
            continue
    return atoms

def _candidate_ionogenic_site_atoms(mol: Chem.Mol) -> list[Chem.Atom]:
    ring_info = mol.GetRingInfo().AtomRings()
    atoms: list[Chem.Atom] = []
    seen: set[int] = set()

    def add(atom: Chem.Atom):
        idx = atom.GetIdx()
        if idx not in seen:
            seen.add(idx)
            atoms.append(atom)
    for atom in _candidate_cation_atoms(mol):
        add(atom)
    for atom in mol.GetAtoms():
        anum = atom.GetAtomicNum()
        if anum not in (7, 15):
            continue
        if atom.GetFormalCharge() > 0:
            add(atom)
            continue
        if anum == 15 and atom.GetTotalDegree() >= 4:
            add(atom)
            continue
        neighbors = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1]
        carbon_neighbors = sum((nb.GetAtomicNum() == 6 for nb in neighbors))
        atom_rings = [ring for ring in ring_info if atom.GetIdx() in ring]
        ring_n_counts = [sum((mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in ring)) for ring in atom_rings]
        aromatic_diaza = atom.GetIsAromatic() and any((len(ring) in (5, 6) and n_count >= 2 for ring, n_count in zip(atom_rings, ring_n_counts)))
        aza_fused = any((len(ring) in (5, 6) and n_count >= 2 for ring, n_count in zip(atom_rings, ring_n_counts)))
        if aromatic_diaza or aza_fused:
            add(atom)
            continue
        if anum == 7 and atom.GetTotalDegree() >= 3 and (carbon_neighbors >= 2):
            add(atom)
            continue
    return atoms


def _aem_for_unit(smiles: str) -> dict[str, float]:
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return {}
    return _mechanism_core_descriptors(mol)


def _public_name(internal_name: str) -> str:
    suffix = internal_name.split("__", 1)[1]
    suffix = suffix.replace("grotthuss_window", "hydration_transport_window")
    suffix = suffix.replace("_dist_", "_distance_")
    if suffix.endswith("_proxy"):
        suffix = suffix[:-6]
    return "aem__" + suffix


def compute_aem_descriptors(
    repeating_unit_1: str,
    repeating_unit_2: str = "",
    ratio_1: float | None = 1.0,
    ratio_2: float | None = 0.0,
) -> dict[str, float]:
    ru1 = canonical_smiles(repeating_unit_1)
    ru2 = canonical_smiles(repeating_unit_2)
    r1, r2, _ = normalize_ratios(ratio_1, ratio_2, ru1, ru2)
    unit1 = _aem_for_unit(ru1)
    unit2 = _aem_for_unit(ru2)
    merged = _composition_aware_merge(unit1, unit2, r1, r2)
    selected = {
        key: value
        for key, value in merged.items()
        if key.startswith(("mech__", "mix__"))
    }
    return {_public_name(key): float(value) for key, value in sorted(selected.items())}
