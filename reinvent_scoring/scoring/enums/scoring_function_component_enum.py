from dataclasses import dataclass


@dataclass(frozen=True)
class ScoringFunctionComponentNameEnum:
    PARALLEL_ROCS_SIMILARITY = "parallel_rocs_similarity"
    SELECTIVITY = "selectivity"
    PREDICTIVE_PROPERTY = "predictive_property"
    CHEMPROP = "chemprop"
    ROCS_SIMILARITY = "rocs_similarity"
    MATCHING_SUBSTRUCTURE = "matching_substructure"
    TANIMOTO_SIMILARITY = "tanimoto_similarity"
    JACCARD_DISTANCE = "jaccard_distance"
    CUSTOM_ALERTS = "custom_alerts"
    QED_SCORE = "qed_score"
    MOLECULAR_WEIGHT = "molecular_weight"
    NUM_ROTATABLE_BONDS = "num_rotatable_bonds"
    NUM_HBD_LIPINSKI = "num_hbd_lipinski"
    NUM_HBA_LIPINSKI = "num_hba_lipinski"
    NUM_RINGS = "num_rings"
    NUM_AROMATIC_RINGS = "num_aromatic_rings"
    NUM_ALIPHATIC_RINGS = "num_aliphatic_rings"
    TPSA = "tpsa"
    SLOGP = "slogp"
    GRAPH_LENGTH = "graph_length"
    NUMBER_OF_STEREO_CENTERS = "number_of_stereo_centers"
    TOTAL_SCORE = "total_score" # there is no actual component corresponding to this type
    REACTION_FILTERS = "reaction_filters"

    # Link invent specific
    LINKER_EFFECTIVE_LENGTH = "linker_effective_length"
    LINKER_GRAPH_LENGTH = "linker_graph_length"
    LINKER_LENGTH_RATIO = "linker_length_ratio"
    LINKER_NUM_RINGS = "linker_num_rings"
    LINKER_NUM_ALIPHATIC_RINGS = "linker_num_aliphatic_rings"
    LINKER_NUM_AROMATIC_RINGS = "linker_num_aromatic_rings"
    LINKER_NUM_SP_ATOMS = "linker_num_sp_atoms"
    LINKER_NUM_SP2_ATOMS = "linker_num_sp2_atoms"
    LINKER_NUM_SP3_ATOMS = "linker_num_sp3_atoms"
    LINKER_NUM_HBA = "linker_num_hba"
    LINKER_NUM_HBD = "linker_num_hbd"
    LINKER_MOL_WEIGHT = "linker_mol_weight"
    LINKER_RATIO_ROTATABLE_BONDS = "linker_ratio_rotatable_bonds"

    #NOTE: components below are AZ specific
    SA_SCORE = "sa_score"
    AZDOCK = "azdock"
    DOCKSTREAM = "dockstream"
    ICOLOS = "icolos"
    AZ_LOGD74_PIP = "azlogd74"
    CACO2_INTR_PIP = "caco2-intrinsic-papp"
    CACO2_EFFLUX_PIP = "caco2-efflux"
    HH_CLINT_PIP = "hh-clint"
    HLM_CLINT_PIP = "hlm-clint"
    RH_CLINT_PIP = "rh-clint"
    SOLUBILITY_DD_PIP = "solubility-dd"
    HERG_PIP = "herg"
    KPUU_PIP = "rat-kpuu-brain"
    RAT_PK_PIP = "rat-pk"
    CLAB_TOP_20 = "clab_top_20"
    RA_SCORE = "rascore"
    AIZYNTH = "aizynth"
    QPTUNA_PIP_MODEL = "optuna-multi"
    THP1_CYTOTOXICITY = "thp1-class"
    GENERAL_REST = "general_rest"

    # Add new RDKit shape components
    RDKIT_SHAPE_SIMILARITY = "rdkit_shape_similarity"
    PARALLEL_RDKIT_SHAPE_SIMILARITY = "parallel_rdkit_shape_similarity"

    # Add Roshambo shape components
    ROSHAMBO_SHAPE_SIMILARITY = "roshambo_shape_similarity"
