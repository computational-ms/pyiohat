"""Ident base parser class."""
import multiprocessing as mp

import ahocorasick
import numpy as np
import pandas as pd
import regex as re
from chemical_composition import ChemicalComposition
from chemical_composition.chemical_composition_kb import PROTON
from loguru import logger
from peptide_mapper.mapper import UPeptideMapper
from unimod_mapper.unimod_mapper import UnimodMapper

from pyprotista.parsers.base_parser import BaseParser
from pyprotista.parsers.misc import (
    get_isotopologue_accuracy,
    init_custom_cc,
)
from pyprotista.utils import merge_and_join_dicts

mp.set_start_method(method="fork")


class IdentBaseParser(BaseParser):
    """Base class of all ident parsers."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.DELIMITER = self.params.get("delimiter", "<|>")
        self.PROTON = PROTON
        self.df = None
        self.mod_mapper = UnimodMapper(xml_file_list=self.xml_file_list)
        self.params["mapped_mods"] = self.mod_mapper.map_mods(
            mod_list=self.params.get("modifications", [])
        )
        successfully_mapped_mods = set(
            mod["name"] for mod in self.params["mapped_mods"]["fix"]
        ) | set(mod["name"] for mod in self.params["mapped_mods"]["opt"])
        self.non_mappable_mods = set(
            mod["name"] for mod in self.params.get("modifications", [])
        ).difference(successfully_mapped_mods)
        self.mod_dict = self._create_mod_dicts()
        self.rt_truncate_precision = 2
        self.reference_dict = {
            "search_engine": None,
            "spectrum_id": None,
            "modifications": None,
            "retention_time_seconds": None,
        }
        self.required_headers = self._load_model("ident_parser_model.json")
        self.col_order = pd.Series(self.required_headers.keys())

    def _calc_mz(self, mass, charge):
        """Calculate mass-to-charge ratio.

        Args:
            mass (pd.Series): masses
            charge (pd.Series): charges
        Returns:
            (pd.Series): m/z
        """
        return (
            mass.astype(float) + (charge.astype(int) * self.PROTON)
        ) / charge.astype(int)

    def _calc_mass(self, mz, charge):
        """Calculate mass from mass-to-charge ratio.

        Args:
            mz (pd.Series): mass-to-charge
            charge (pd.Series): charges
        Returns:
            (pd.Series): mass
        """
        return mz.astype(float) * charge.astype(int) - (
            charge.astype(int) * self.PROTON
        )

    def _create_mod_dicts(self):
        """
        Create dict containing meta information about static and variable mods.

        Returns:
            mod_dict (dict): mapped modifications and information
        """
        mod_dict = {}
        for mod_type in ["fix", "opt"]:
            for modification in self.params["mapped_mods"][mod_type]:
                aa = modification["aa"]
                pos = modification["position"]
                name = modification["name"]
                if name not in mod_dict.keys():
                    mod_dict[name] = {
                        "mass": modification["mass"],
                        "aa": set(),
                        "position": set(),
                    }
                mod_dict[name]["aa"].add(aa)

                mod_dict[name]["aa"].add(pos)
                mod_dict[name]["position"].add(pos)

        return mod_dict

    def clean_up_modifications(self):
        """Sanitizes modstrings generated by engine parsers.

        Modifications are sorted by position and leading, repeated or trailing delimiters are removed
        Operations are performed inplace on self.df
        """
        # Remove any trailing or leading delimiters or only-delimiter modstrings
        self.df.loc[:, "modifications"] = self.df.loc[:, "modifications"].str.replace(
            r"^;+(?=\w)", "", regex=True
        )
        self.df.loc[:, "modifications"] = self.df.loc[:, "modifications"].str.replace(
            r"(?<=\w);+$", "", regex=True
        )
        self.df.loc[:, "modifications"] = self.df.loc[:, "modifications"].str.replace(
            r"^;+$", "", regex=True
        )

        # Ensure same order of modifications
        self.df.loc[:, "modifications"] = (
            self.df["modifications"].fillna("").str.split(";").apply(self.sort_mods)
        )

    def sort_mods(self, data):
        """Sort mods based on position in peptide and ensures unique mod:pos pairs.

        Args:
            data (list): list of modifications with style "Mod:position"
        Returns:
            sorted_formatted_mods (str): String with sorted mods in style "Mod1:pos1;Modn:posn"
        """
        data = set(data)
        sort_pattern = r"([\w\-\(\)\>\:]+)(?:\:)(\d+)"
        positions = [int(re.search(sort_pattern, d).group(2)) for d in data if d != ""]
        names = [re.search(sort_pattern, d).group(1) for d in data if d != ""]
        sorted_mods = sorted(zip(names, positions), key=lambda x: x[1])
        sorted_mods_str = [[str(i) for i in m] for m in sorted_mods]
        sorted_formatted_mods = ";".join([":".join(m) for m in sorted_mods_str])
        return sorted_formatted_mods

    def assert_only_iupac_and_missing_aas(self):
        """Assert that only IUPAC nomenclature one letter amino acids are used in sequence.

        Non-IUPAC designations are dropped (except for selenocysteine).
        Operations are performed inplace.
        """
        self.df["sequence"] = self.df["sequence"].str.upper()
        # Added X for missing AAs
        iupac_aas = set("ACDEFGHIKLMNPQRSTUVWY")
        iupac_conform_seqs = self.df["sequence"].apply(
            lambda seq: set(seq).issubset(iupac_aas)
        )
        if any(~iupac_conform_seqs):
            self.df = self.df.loc[iupac_conform_seqs, :]
            logger.warning(
                f"sequences are not IUPAC conform. {(~iupac_conform_seqs).sum()} PSMs were dropped."
            )

    def add_protein_ids(self):
        """Add all Protein IDs that matching the sequence.

        Operations are performed inplace on self.df
        """
        peptide_mapper = UPeptideMapper(self.params["database"])
        mapped_peptides = peptide_mapper.map_peptides(self.df["sequence"].tolist())

        peptide_mappings = [
            merge_and_join_dicts(mapped_peptides[seq], self.DELIMITER)
            for seq in self.df["sequence"]
        ]

        columns_translations = {
            "start": "sequence_start",
            "end": "sequence_stop",
            "post": "sequence_post_aa",
            "id": "protein_id",
            "pre": "sequence_pre_aa",
        }
        new_columns = pd.DataFrame(peptide_mappings)
        new_columns.rename(columns=columns_translations, inplace=True)

        self.df.loc[:, new_columns.columns] = new_columns.values
        new_columns = new_columns.dropna(axis=0, how="all")
        if len(new_columns) != len(self.df):
            logger.warning(
                f"{len(self.df)-len(new_columns)} PSMs were dropped because their respective sequences could not be mapped."
            )
        self.df = self.df.iloc[new_columns.index, :].reset_index(drop=True)

    def check_enzyme_specificity(self):
        """Check consistency of N/C-terminal cleavage sites.

        Calculates number of missed cleavage sites.
        Operations are performed inplace.
        """
        if self.params["enzyme"] == ".^":
            self.df.loc[:, ["enzn", "enzc"]] = True
            self.df.loc[:, "missed_cleavages"] = 0
            return None

        enzyme_pattern = self.params["enzyme"]
        integrity_strictness = self.params["terminal_cleavage_site_integrity"]

        pren_seq = (
            pd.concat(
                [
                    self.df["sequence_pre_aa"].str.split(rf"{self.DELIMITER}"),
                    self.df["sequence"].str[:1],
                ],
                axis=1,
            )
            .explode("sequence_pre_aa")
            .sum(axis=1)
        )
        self.df.loc[:, "enzn"] = (
            pren_seq.str.split(rf"{enzyme_pattern}").str[0].str.len() == 1
        ).groupby(pren_seq.index).agg(integrity_strictness) | (
            pren_seq.str[0] == "-"
        ).groupby(
            pren_seq.index
        ).agg(
            integrity_strictness
        )
        postc_seq = (
            pd.concat(
                [
                    self.df["sequence"].str[-1:],
                    self.df["sequence_post_aa"].str.split("<\\|>"),
                ],
                axis=1,
            )
            .explode("sequence_post_aa")
            .sum(axis=1)
        )
        self.df.loc[:, "enzc"] = (
            postc_seq.str.split(rf"{enzyme_pattern}").str[0].str.len() == 1
        ).groupby(postc_seq.index).agg(integrity_strictness) | (
            postc_seq.str[-1] == "-"
        ).groupby(
            postc_seq.index
        ).agg(
            integrity_strictness
        )

        internal_cuts = self.df["sequence"].str.split(rf"{enzyme_pattern}")
        self.df.loc[:, "missed_cleavages"] = (
            internal_cuts.apply(len)
            - internal_cuts.apply(lambda row: "" in row).astype(int)
            - 1
        )

    def calc_masses_offsets_and_composition(self):
        """Calculate chemical composition theoretical masses, and mass-to-charge ratios.

        Offsets are calculated between theoretical and experimental mass-to-charge ratio.
        Operations are performed inplace on self.df
        """
        all_compositions = {}
        iupac_aas = tuple("ACDEFGHIKLMNPQRSTUVWY")
        cc = ChemicalComposition(
            unimod_file_list=self.params.get("xml_file_list", None)
        )
        for aa in iupac_aas:
            cc.use(sequence=aa)
            all_compositions[aa] = cc.copy()
        for mod in self.mod_dict.keys():
            all_compositions[mod] = self.mod_mapper.name_to_composition(mod)[0]
        elements = list(
            set(element for comp in all_compositions.values() for element in comp)
        )
        # C and H first then the rest of the elements
        elements = list(
            sorted(
                elements, key=lambda e: "00" if e == "C" else "01" if e == "H" else e
            )
        )
        atom_counts = np.zeros(shape=(len(self.df), len(elements)), dtype=int)
        sequences = self.df["sequence"].to_numpy(dtype=str)
        modifications = self.df["modifications"].to_numpy(dtype=str)
        for aa in iupac_aas:
            ordered_element_multiplier = []
            for element in elements:
                ordered_element_multiplier += [all_compositions[aa].get(element, 0)]
            ordered_element_multiplier = np.array(ordered_element_multiplier)
            atom_counts += np.outer(
                np.char.count(sequences, aa), ordered_element_multiplier
            )
        # Remove water (peptide bonds)
        water = np.zeros(shape=(1, len(elements)), dtype=int)
        water[0, elements.index("H")] = 2
        water[0, elements.index("O")] = 1
        atom_counts -= np.outer(np.maximum(np.char.str_len(sequences) - 1, 0), water)
        for mod in self.mod_dict.keys():
            ordered_element_multiplier = []
            for element in elements:
                ordered_element_multiplier += [all_compositions[mod].get(element, 0)]
            ordered_element_multiplier = np.array(ordered_element_multiplier)
            atom_counts += np.outer(
                np.char.count(modifications, mod + ":"), ordered_element_multiplier
            )
        for i, element in enumerate(elements):
            if i == 0:
                chemical_compositions = np.char.add(
                    elements[i] + "(", atom_counts.astype(str)[:, i]
                )
            else:
                next_element_string = np.char.add(
                    ")" + elements[i] + "(", atom_counts.astype(str)[:, i]
                )
                chemical_compositions = np.char.add(
                    chemical_compositions, next_element_string
                )
            if i == len(elements) - 1:
                chemical_compositions = np.char.add(chemical_compositions, ")")
        self.df["chemical_composition"] = chemical_compositions
        # Remove empties
        self.df["chemical_composition"] = self.df["chemical_composition"].str.replace(
            r"\d*[A-Z][a-z]?\(0\)", "", regex=True
        )

        isotope_mass_lookup = {}
        for element, isotope_data in cc.isotopic_distributions.items():
            if isinstance(isotope_data, dict):
                # this has extra data
                isotope_list = []
                for fraction in isotope_data:
                    isotope_list.extend(isotope_data[fraction])
            else:
                isotope_list = isotope_data

            for iso_abund in isotope_list:
                isotope_mass_key = f"{round(iso_abund[0])}{element}"
                isotope_mass_lookup[isotope_mass_key] = iso_abund[0]

        monoisotopic_element_masses = []
        for element in elements:
            try:
                mass = max(cc.isotopic_distributions[element], key=lambda d: d[1])[0]
            except KeyError:
                # Element is an isotope
                mass = isotope_mass_lookup[element]
            monoisotopic_element_masses.extend([mass])
        monoisotopic_element_masses = np.array(monoisotopic_element_masses)
        self.df["ucalc_mass"] = (atom_counts * monoisotopic_element_masses).sum(axis=1)

        with mp.Pool(
            self.params.get("cpus", mp.cpu_count() - 1),
            initializer=init_custom_cc,
            initargs=(
                get_isotopologue_accuracy,
                self.PROTON,
            ),
        ) as pool:
            acc = pool.starmap(
                get_isotopologue_accuracy,
                zip(
                    self.df["chemical_composition"].values,
                    self.df["charge"].astype(int).values,
                    self.df["exp_mz"].values,
                ),
                chunksize=1,
            )
        self.df.loc[:, "accuracy_ppm"] = acc
        self.df.loc[:, "ucalc_mz"] = self._calc_mz(
            mass=self.df["ucalc_mass"], charge=self.df["charge"]
        )
        self.df.loc[:, "accuracy_ppm_C12"] = (
            (self.df["exp_mz"].astype(float) - self.df["ucalc_mz"])
            / self.df["ucalc_mz"]
            * 1e6
        )
        # Clear rows with non-mappable mods
        for mod in self.non_mappable_mods:
            self.df.loc[
                self.df["modifications"].str.contains(mod),
                (
                    "chemical_composition",
                    "ucalc_mass",
                    "ucalc_mz",
                    "accuracy_ppm",
                    "accuracy_ppm_C12",
                ),
            ] = np.nan

    def _read_meta_info_lookup_file(self):
        """Read meta info lookup file.

        Returns:
            rt_lookup (dict of int: dict): dict with spectrum ids as top level key
                                           values are new dicts with all rt values as keys
                                           values are lists with [file, precursor_mz]
        """
        rt_lookup = pd.read_csv(self.params["rt_pickle_name"], compression="infer")
        rt_lookup["rt_unit"] = rt_lookup["rt_unit"].replace(
            {"second": 1, "minute": 60, "s": 1, "min": 60}
        )
        rt_lookup["retention_time_seconds"] = (
            rt_lookup[["rt", "rt_unit"]].astype(float).product(axis=1)
        )
        rt_lookup.drop(columns=["rt", "rt_unit"], inplace=True)

        rt_lookup.set_index(
            [
                "spectrum_id",
                "retention_time_seconds",
            ],
            inplace=True,
        )
        rt_lookup = rt_lookup.loc[:, ["file", "precursor_mz"]]
        rt_lookup = (
            rt_lookup.groupby(level=0)
            .apply(lambda grp: grp.xs(grp.name).T.to_dict("list"))
            .to_dict()
        )
        return rt_lookup

    def get_meta_info(self):
        """Extract meta information.

        Experimental mass-to-charge ratios, retention times, file names,
        and spectrum titles are added.
        Operations are performed inplace on self.df
        """
        rt_lookup = self._read_meta_info_lookup_file()
        self.df["spectrum_id"] = self.df["spectrum_id"].astype(int)
        if self.style in ("comet_style_1", "omssa_style_1"):
            logger.warning(
                "This engine does not provide retention time information. Grouping only by Spectrum ID. This may cause problems when working with multi-file inputs."
            )
            for name, grp in self.df.groupby("spectrum_id"):
                mappable_within_precision = list(rt_lookup[name].keys())
                if len(mappable_within_precision) == 1:
                    self.df.loc[
                        grp.index,
                        ("raw_data_location", "exp_mz", "retention_time_seconds"),
                    ] = rt_lookup[name][mappable_within_precision[0]] + [
                        mappable_within_precision[0]
                    ]
                else:
                    logger.error(
                        f"Could not uniquely assign meta data to spectrum id {name}."
                    )
        else:
            self.df["retention_time_seconds"] = self.df[
                "retention_time_seconds"
            ].astype(float)
            for name, grp in self.df.groupby(["spectrum_id", "retention_time_seconds"]):
                meta_rts = rt_lookup[name[0]].keys()
                mappable_within_precision = [
                    rt
                    for rt in meta_rts
                    if abs(name[1] - rt) <= 10**-self.rt_truncate_precision
                ]
                if len(mappable_within_precision) == 1:
                    self.df.loc[
                        grp.index,
                        ("raw_data_location", "exp_mz", "retention_time_seconds"),
                    ] = rt_lookup[name[0]][mappable_within_precision[0]] + [
                        mappable_within_precision[0]
                    ]
                else:
                    logger.error(
                        f"Could not uniquely assign meta data to spectrum id, retention time {name}."
                    )

        self.df.loc[:, "spectrum_title"] = (
            self.df["raw_data_location"]
            + "."
            + self.df["spectrum_id"].astype(str)
            + "."
            + self.df["spectrum_id"].astype(str)
            + "."
            + self.df["charge"].astype(str)
        )

    def add_ranks(self):
        """Ranks are calculated based on the engine scoring column at Spectrum ID level.

        Operations are performed inplace on self.df
        """
        eng_name = self.df["search_engine"].unique()[0]
        score_col = self.params["validation_score_field"][eng_name]
        top_is_highest = self.params["bigger_scores_better"][eng_name]
        ranking_needs_to_be_ascending = False if top_is_highest is True else True
        self.df.loc[:, score_col] = self.df[score_col].astype(float)
        self.df.loc[:, "rank"] = (
            self.df.groupby("spectrum_id")[score_col]
            .rank(ascending=ranking_needs_to_be_ascending, method="min")
            .astype(int)
        )

    def add_decoy_identity(self):
        """Add boolean decoy state if designated decoy prefix is in Protein IDs.

        Also marks peptides which were assigned decoy but were immutable during
        target decoy generation.
        Operations are performed inplace on self.df
        """
        decoy_tag = self.params.get("decoy_tag", "decoy_")
        self.df.loc[:, "is_decoy"] = self.df["protein_id"].str.contains(decoy_tag)
        if self.immutable_peptides is not None:
            auto = ahocorasick.Automaton()
            for seq in self.immutable_peptides:
                auto.add_word(seq, seq)
            auto.make_automaton()
            self.df.loc[:, "is_immutable"] = [
                sum([len(match) for _, match in auto.iter_long(seq)]) == len(seq)
                for seq in self.df["sequence"]
            ]
        else:
            self.df.loc[:, "is_immutable"] = False

    def process_unify_style(self):
        """Combine all additional operations that are needed to calculate new columns and sanitize the dataframe.

        Operations are performed inplace on self.df
        """
        self.df.drop_duplicates(inplace=True, ignore_index=True)
        self.clean_up_modifications()
        self.assert_only_iupac_and_missing_aas()
        self.add_protein_ids()
        self.get_meta_info()
        self.calc_masses_offsets_and_composition()
        self.check_enzyme_specificity()
        self.add_ranks()
        self.add_decoy_identity()
        self.sanitize()
