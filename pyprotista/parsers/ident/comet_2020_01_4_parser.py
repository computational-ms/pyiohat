"""Engine parser."""

import numpy as np
import pandas as pd
import regex as re
import xml.etree.cElementTree as etree

from pyprotista.parsers.ident_base_parser import IdentBaseParser

element_tag_prefix = "{http://psidev.info/psi/pi/mzIdentML/1.2}"


def _iterator_xml(xml_file, mapping_dict):
    version = ""
    stage = 0

    modifications = {}
    sequence = {}
    peptide_lookup = {}

    spec_results = {}
    spec_ident_items = []
    spec_records = []

    modification_mass_map = {}
    mod_name = ""
    fixed_mods = {}

    for event, entry in etree.iterparse(xml_file):
        entry_tag = entry.tag

        if entry_tag == (f"{element_tag_prefix}PeptideEvidence"):
            stage = 0
            continue
            # Back to 0, so no cvParam gets written
        elif entry_tag == (f"{element_tag_prefix}ModificationParams"):
            stage = 0
        elif entry_tag == (f"{element_tag_prefix}PeptideEvidenceRef"):
            continue

        elif stage == 1:
            sequence, modifications, peptide_lookup = _peptide_lookup(
                entry, entry_tag, sequence, modifications, peptide_lookup
            )
        elif stage == 2:
            modification_mass_map, mod_name, fixed_mods = _get_modifications(
                entry, entry_tag, modification_mass_map, mod_name, fixed_mods
            )
        elif stage == 3:
            spec_results, spec_ident_items, spec_records = _spec_records(
                entry,
                entry_tag,
                spec_results,
                spec_ident_items,
                spec_records,
                mapping_dict,
            )

        elif entry_tag == (f"{element_tag_prefix}DBSequence"):
            stage = 1
            continue
            # Short before peptide_lookup begins
        elif entry_tag == (f"{element_tag_prefix}AdditionalSearchParams"):
            stage = 2
            continue
        elif entry_tag == (f"{element_tag_prefix}Inputs"):
            stage = 3
            continue
            # Short before spec_record begins

        elif entry_tag == (f"{element_tag_prefix}AnalysisSoftware"):
            version = "comet_" + "_".join(
                re.findall(r"([/d]*\d+)", entry.attrib["version"])
            )

        entry.clear()
        if entry_tag == (f"{element_tag_prefix}SpectrumIdentificationList"):
            return (
                version,
                peptide_lookup,
                spec_records,
                modification_mass_map,
                fixed_mods,
            )


def _peptide_lookup(entry, entry_tag, sequence, modifications, peptide_lookup):
    if entry_tag == (f"{element_tag_prefix}PeptideSequence"):
        sequence = {"sequence": entry.text}
    elif entry_tag == (f"{element_tag_prefix}Modification"):
        modifications.update(
            {"monoisotopicMassDelta": entry.attrib["monoisotopicMassDelta"]}
        )
        modifications.update({"location": entry.attrib["location"]})
    elif entry_tag == (f"{element_tag_prefix}Peptide"):
        if len(modifications) > 0:
            peptide_lookup[entry.attrib["id"]] = {"modifications": modifications}
        else:
            peptide_lookup[entry.attrib["id"]] = {"modifications": ""}
        peptide_lookup[entry.attrib["id"]].update(sequence)
        modifications = {}
        sequence = ""
    return sequence, modifications, peptide_lookup


def _get_modifications(entry, entry_tag, modification_mass_map, mod_name, fixed_mods):
    if entry_tag == (f"{element_tag_prefix}cvParam"):
        mod_name = entry.attrib["name"]
    elif entry_tag == (f"{element_tag_prefix}SearchModification"):
        modification_mass_map[entry.attrib["massDelta"]] = mod_name
        if entry.attrib["fixedMod"] == "true":
            fixed_mods.update({entry.attrib["residues"]: mod_name})
    return modification_mass_map, mod_name, fixed_mods


def _spec_records(
    entry, entry_tag, spec_results, spec_ident_items, spec_records, mapping_dict
):
    if entry_tag == (f"{element_tag_prefix}cvParam"):
        if entry.attrib["name"] in mapping_dict:
            spec_results.update(
                {mapping_dict[entry.attrib["name"]]: entry.attrib["value"]}
            )
    elif entry_tag == (f"{element_tag_prefix}SpectrumIdentificationItem"):
        for i in list(entry.attrib):
            if i in mapping_dict:
                spec_results.update({mapping_dict[i]: entry.attrib[i]})
        spec_ident_items.append(spec_results)
        spec_results = {}
    elif entry_tag == (f"{element_tag_prefix}SpectrumIdentificationResult"):
        for i in spec_ident_items:
            i.update(spec_results)
            i.update({"spectrum_id": entry.attrib["spectrumID"].lstrip("scan=")})
            spec_records.append(i)
        spec_results = {}
        spec_ident_items = []
    return spec_results, spec_ident_items, spec_records


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
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})

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

    def _map_mods_and_sequences(
        self, fixed_mods, modification_mass_map, peptide_lookup
    ):
        if len(fixed_mods) > 0:
            fixed_mod_strings = []
            for fm_res, fm_name in fixed_mods.items():
                fixed_mod_strings.append(
                    self.df["sequence"]
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
                pd.concat(fixed_mod_strings, axis=1)
                .agg(";".join, axis=1)
                .str.rstrip(";")
            )
        for i in peptide_lookup.keys():
            if len(peptide_lookup[i]["modifications"]) != 0:
                monoisotopicMassDelta = peptide_lookup[i]["modifications"][
                    "monoisotopicMassDelta"
                ]
                location = peptide_lookup[i]["modifications"]["location"]
                peptide_lookup[i][
                    "modifications"
                ] = f"{modification_mass_map[monoisotopicMassDelta]}:{location}"

        seq_mods = pd.DataFrame(self.df["sequence"].map(peptide_lookup).to_list())
        self.df.loc[:, "modifications"] = (
            seq_mods["modifications"].str.cat(fixed_mod_strings, sep=";").str.strip(";")
        )
        self.df.loc[:, "sequence"] = seq_mods["sequence"]

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        (
            version,
            peptide_lookup,
            spec_records,
            modification_mass_map,
            fixed_mods,
        ) = _iterator_xml(self.input_file, self.mapping_dict)
        self.df = pd.DataFrame(spec_records)
        self.df["modifications"] = None
        self.df["retention_time_seconds"] = None
        self.df["search_engine"] = version
        self._map_mods_and_sequences(fixed_mods, modification_mass_map, peptide_lookup)
        self.process_unify_style()

        return self.df
