"""Quant parser."""

import re
from pathlib import Path

import pandas as pd
from chemical_composition import ChemicalComposition

from pyiohat.parsers.quant_base_parser import QuantBaseParser
from pyiohat.parsers.misc import get_compositions_and_monoisotopic_masses


class FlashLFQ_1_2_0_Parser(QuantBaseParser):
    """File parser for Flash LFQ."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "flash_lfq_style_1"
        self.mapping_dict = {
            "Peak Apex Mass Error (ppm)": "accuracy_mz",
            "Precursor Charge": "charge",
            "Base Sequences Mapped": "flashlfq:base_sequences_mapped",
            "Full Sequence": "flashlfq:full_sequence",
            "Full Sequences Mapped": "flashlfq:full_sequences_mapped",
            "MBR Score": "flashlfq:mbr_score",
            "MS2 Retention Time": "flashlfq:ms2_retention_time",
            "Num Charge States Observed": "flashlfq:num charge states observed",
            "Peak Charge": "flashlfq:peak_charge",
            "Peak Detection Type": "flashlfq:peak_detection_type",
            "Peak RT Apex": "flashlfq:peak_rt_apex",
            "Peak RT End": "flashlfq:peak_rt_end",
            "Peak RT Start": "flashlfq:peak_rt_start",
            "Peak Split Valley RT": "flashlfq:peak_split_valley_rt",
            "Peptide Monoisotopic Mass": "flashlfq:peptide_monoisotopic_mass",
            "Protein Group": "flashlfq:protein_group",
            "PSMs Mapped": "flashlfq:psms_mapped",
            "Peak intensity": "quant_value",
            "File Name": "raw_data_location",
            "Peak MZ": "reported_mz",
            "Theoretical MZ": "theoretical_mz",
            "Base Sequence": "trivial_name",
        }
        self.df = pd.read_csv(self.input_file, delimiter="\t")
        self.df.rename(columns=self.mapping_dict, inplace=True)
        self.round_precision = 5

        self.rt_to_spec_id = {}
        self.filestem_to_path = {}
        for spec_id, metadata in self.rt_lookup.items():

            file_path = metadata.get("lineage_root")
            rt_value = metadata.get("rt")
            if isinstance(file_path, list):
                file_path = file_path[0] if file_path else None
            if file_path is not None and rt_value is not None:
                rt_list = rt_value if isinstance(rt_value, list) else [rt_value]
                for rt in rt_list:
                    rounded_rt = round(rt, self.round_precision)
                self.rt_to_spec_id.setdefault(file_path, {})[rounded_rt] = spec_id
                stem = Path(file_path).stem
                if stem not in self.filestem_to_path:
                    self.filestem_to_path[stem] = file_path

            # round(list(d.keys())[0], self.round_precision): key
            # for key, d in self.rt_lookup.items()
        # }
        self.cc = ChemicalComposition()
        self.IUPAC_AAS = tuple("ACDEFGHIKLMNPQRSTUVWY")
        self.metadata = {
            "File Origin": "Flashlfq",
            "Parser": "pyiohat/parsers/quant/flash_lfq_1_2_0_parser.py",
        }

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_tsv = file.as_posix().endswith(".tsv")
        flash_lfq_columns = {
            "File Name",
            "Base Sequence",
            "Full Sequence",
            "Protein Group",
            "Peptide Monoisotopic Mass",
            "MS2 Retention Time",
            "Precursor Charge",
            "Theoretical MZ",
            "Peak intensity",
            "Peak RT Start",
            "Peak RT Apex",
            "Peak RT End",
            "Peak MZ",
            "Peak Charge",
            "Num Charge States Observed",
            "Peak Detection Type",
            "MBR Score",
            "PSMs Mapped",
            "Base Sequences Mapped",
            "Full Sequences Mapped",
            "Peak Split Valley RT",
            "Peak Apex Mass Error (ppm)",
        }
        with open(file.as_posix()) as f:
            head = set(f.readline().replace("\n", "").split("\t"))
        headers_match = len(flash_lfq_columns.difference(head)) == 0
        return is_tsv and headers_match

    def unify(self):
        """Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        # TODO: fill this
        # raise NotImplementedError
        # do column conversion here
        self.df["spectrum_id"] = -1
        self.df["ident_reference"] = -1
        self.df["retention_time_seconds"] = -1
        self.df["accuracy_ppm_C12"] = -1
        self.df["quant_score"] = -1
        self.df["fwhm"] = -1
        self.df["label"] = "LabelFree"
        self.df["quant_group"] = ""
        self.df["ccs"] = -1
        self.df.loc[self.df["flashlfq:mbr_score"].isna(), "flashlfq:mbr_score"] = -1

        self.get_chemical_composition()
        self.get_meta_info()
        self.calculate_accuracy()
        self.df["quant_run_id"] = "FlashLFQ"
        self.df["quant_engines"] = "FlashLFQ_1_2_0"

        self.process_unify_style()
        return self.df

    def get_meta_info(self):
        rounded_rts = (self.df["flashlfq:ms2_retention_time"]).apply(
            round, args=(self.round_precision,)
        )
        rounded_rts = pd.DataFrame(
            {"file": self.df["raw_data_location"], "rt": rounded_rts}
        )

        self.df["ident_reference"] = [
            self.rt_to_spec_id[f][r] for i, (f, r) in rounded_rts.iterrows()
        ]

    def get_chemical_composition(self):
        mods = self.df["flashlfq:full_sequence"].apply(self.translate_mods)
        seqs = self.df["trivial_name"]

        # move to base parser?
        all_compositions = {}
        for aa in self.IUPAC_AAS:
            self.cc.use(sequence=aa)
            all_compositions[aa] = self.cc.copy()
        for mod in self.mapped_mod_names:
            all_compositions[mod] = self.mod_mapper.name_to_composition(mod)[0]

        compositions, mono_masses = get_compositions_and_monoisotopic_masses(
            sequences=seqs.to_numpy(dtype=str),
            modifications=mods.to_numpy(dtype=str),
            compositions=all_compositions,
            isotopic_distributions=self.cc.isotopic_distributions,
        )
        self.df.loc[:, "chemical_composition"] = compositions

    def translate_mods(self, full_sequence):
        """Extract modifications from full_sequence and format as {mod_1}:{pos_1};{mod_n}:{pos_n}.

        Args:
            full_sequence (str): sequence including mod (e.g. ELVISC[Carbamidomethyl]M[Oxidation])

        Returns:
            str: extracted mods as described in summary
        """
        # TODO extract C-terminal mods, position should be seq_len + 1
        cumulative_match_length = 0
        regex = re.compile(r"\[(.*?)\]")
        mods = []
        for match in regex.finditer(full_sequence):
            mods.append(f"{match.group(1)}:{match.start() - cumulative_match_length}")
            cumulative_match_length += len(match.group())
        return ";".join(mods)
