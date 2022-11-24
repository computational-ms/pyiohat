"""Engine parser."""
import numpy as np
import pandas as pd
import regex as re
import xml.etree.ElementTree as etree
from loguru import logger
from pyprotista.parsers.ident_base_parser import IdentBaseParser


def get_modification_mass_mapping(file):
    """Iterate over file to get information on modifications.

    Args:
        file (str): path to input file

    Returns:
        fixed_mods (dict): residues and corresponding names of modifications, that are fixed
        modification_mass_map (dict): mapping of masses to modification names

    """
    fixed_mods = {}
    modification_mass_map = {}

    # Temporary variables, get overwritten multiple times during iteration.
    stage = 0
    mod_name = ""

    for event, entry in etree.iterparse(file):
        entry_tag = entry.tag
        if entry_tag.endswith("AdditionalSearchParams"):
            stage = 1
        elif entry_tag.endswith("ModificationParams"):
            return fixed_mods, modification_mass_map
        elif stage == 1:
            if entry_tag.endswith("cvParam"):
                mod_name = entry.attrib["name"]
            elif entry_tag.endswith("SearchModification"):
                modification_mass_map[entry.attrib["massDelta"]] = mod_name
                if entry.attrib["fixedMod"] == "true":
                    fixed_mods.update({entry.attrib["residues"]: mod_name})
        entry.clear()


def get_xml_data(file, mapping_dict, fixed_mods, modification_mass_map):
    """Iterate over file to get information on specs.

    Args:
        file (str): path to input file
        mapping_dict (dict): mapping of engine level column names to pyprotista unified column names
        fixed_mods (dict): residues and corresponding names of modifications, that are fixed
        modification_mass_map (dict): mapping of masses to modification names

    Returns:
        version (str): file version
        spec_records (list): spectrum information
    """
    version = ""
    spec_records = []

    # Temporary variables, get overwritten multiple times during iteration
    stage = 0
    modifications = ""
    sequence = {}
    spec_results = {}
    spec_ident_items = []
    peptide_lookup = {}

    for event, entry in etree.iterparse(file):
        entry_tag = entry.tag
        if entry_tag.endswith("PeptideEvidence"):
            # Back to 0, so no cvParam gets written
            stage = 0
        elif stage == 1:
            sequence, modifications, peptide_lookup = get_peptide_lookup(
                entry=entry,
                entry_tag=entry_tag,
                sequence=sequence,
                modifications=modifications,
                peptide_lookup=peptide_lookup,
                fixed_mods=fixed_mods,
                modification_mass_map=modification_mass_map,
            )
        elif stage == 2:
            spec_results, spec_ident_items, spec_records = get_spec_records(
                entry=entry,
                entry_tag=entry_tag,
                spec_results=spec_results,
                spec_ident_items=spec_ident_items,
                spec_records=spec_records,
                mapping_dict=mapping_dict,
                peptide_lookup=peptide_lookup,
            )
        elif entry_tag.endswith("DBSequence"):
            stage = 1
        elif entry_tag.endswith("Inputs"):
            stage = 2
        elif entry_tag.endswith("AnalysisSoftware"):
            version = "comet_" + "_".join(re.findall("[0-9]+", entry.attrib["version"]))
        elif entry_tag.endswith("cvList"):
            if entry_tag != "{http://psidev.info/psi/pi/mzIdentML/1.2}cvList":
                logger.warning("Wrong mzIdentML format - Parser might not operate correctly!")
        entry.clear()
        if entry_tag.endswith("SpectrumIdentificationList"):
            return (
                version,
                spec_records,
            )


def get_peptide_lookup(
    entry,
    entry_tag,
    sequence,
    modifications,
    peptide_lookup,
    fixed_mods,
    modification_mass_map,
):
    """Take one entry at a time to get peptide ids with their corresponding sequences and modifications.

    Args:
        entry (element): xml element
        entry_tag (str): xml tag
        sequence (str): current sequence
        modifications (str): current modifications
        peptide_lookup (dict): peptide id with corresponding modifications and sequence
        fixed_mods (dict): residues and corresponding names of modifications, that are fixed
        modification_mass_map (dict): mapping of masses to modification names

    Returns:
        sequence (str): current sequence
        modifications (str): current modifications
        peptide_lookup (dict): peptide id with corresponding modifications and sequence
    """
    if entry_tag.endswith("PeptideSequence"):
        sequence = entry.text
        if len(fixed_mods) > 0:
            sequence_mod_map = map_mods_sequences(fixed_mods, sequence)
            modifications = sequence_mod_map + ";"
        else:
            modifications = ""
    elif entry_tag.endswith("Modification"):
        mass = entry.attrib["monoisotopicMassDelta"]
        location = entry.attrib["location"]
        mass_name = modification_mass_map[mass]
        modifications += mass_name + ":" + location + ";"
    elif entry_tag.endswith("Peptide"):
        modifications = modifications.rstrip(";").lstrip(";")
        peptide_lookup[entry.attrib["id"]] = (modifications, sequence)
    return sequence, modifications, peptide_lookup


def get_spec_records(
    entry,
    entry_tag,
    spec_results,
    spec_ident_items,
    spec_records,
    mapping_dict,
    peptide_lookup,
):
    """Take one entry at a time to get single specs.

    Args:
        entry (element) : xml element
        entry_tag (str): xml tag
        spec_results (dict): single PSM
        spec_ident_items (list): potentially multiple PSMs
        spec_records (list): information on PSMs
        mapping_dict (dict): mapping of engine level column names to pyprotista unified column names
        peptide_lookup (dict): peptide id with corresponding modifications and sequence

    Returns:
        spec_results (dict): single PSM
        spec_ident_items (list): potentially multiple PSMs
        spec_records (list): information on PSMs
    """
    if entry_tag.endswith("cvParam"):
        if entry.attrib["name"] in mapping_dict:
            spec_results.update(
                {mapping_dict[entry.attrib["name"]]: entry.attrib["value"]}
            )
    elif entry_tag.endswith("SpectrumIdentificationItem"):
        for attribute in list(entry.attrib):
            if attribute in mapping_dict.keys():
                spec_results.update({mapping_dict[attribute]: entry.attrib[attribute]})
        mods, sequence = peptide_lookup[spec_results["sequence"]]
        spec_results.update({"modifications": mods, "sequence": sequence})
        spec_ident_items.append(spec_results)
        spec_results = {}
    elif entry_tag.endswith("SpectrumIdentificationResult"):
        for spec_item in spec_ident_items:
            spec_item.update(spec_results)
            spec_item.update(
                {"spectrum_id": entry.attrib["spectrumID"].lstrip("scan=")}
            )
            spec_records.append(spec_item)
        spec_results = {}
        spec_ident_items = []
    return spec_results, spec_ident_items, spec_records


def map_mods_sequences(fixed_mods, sequence):
    """Map fixed_mods and sequence to get corresponding modifications.

    Args:
        fixed_mods (dict): residues and corresponding names of modifications, that are fixed
        sequence (str): current sequence
    Returns:
        fixed_mod_strings (str): modifications of corresponding sequence
    """
    if len(fixed_mods) == 0:
        logger.error("No fixed mods to map. Shouldn't be here!")
    fixed_mod_strings = []
    for fm_res, fm_name in fixed_mods.items():
        fixed_mod_strings.append(
            pd.Series(sequence)
            .str.split(fm_res)
            .apply(
                lambda l: ";".join(
                    [
                        fm_name + ":" + ind
                        for ind in (
                            np.cumsum(list(map(len, l[:-1]))) + range(1, len(l))
                        ).astype(str)
                    ]
                )
            )
        )

    fixed_mod_strings = (
        pd.concat(fixed_mod_strings, axis=1).agg(";".join, axis=1).str.rstrip(";")
    )
    fixed_mod_strings = fixed_mod_strings[0]
    return fixed_mod_strings


class Comet_2020_01_4_Parser(IdentBaseParser):
    """File parser for Comet."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "comet_style_1"
        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
        }

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_mzid = file.as_posix().endswith(".mzid")

        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for _ in range(10)])
            except StopIteration:
                head = ""
        contains_engine = "Comet" in head
        return is_mzid and contains_engine

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        (fixed_mods, modification_mass_map) = get_modification_mass_mapping(
            self.input_file
        )
        (version, spec_records,) = get_xml_data(
            file=self.input_file,
            mapping_dict=self.mapping_dict,
            fixed_mods=fixed_mods,
            modification_mass_map=modification_mass_map,
        )
        self.df = pd.DataFrame(spec_records)
        self.df["search_engine"] = version
        self.process_unify_style()

        return self.df
