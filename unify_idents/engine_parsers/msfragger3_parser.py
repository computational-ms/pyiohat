from unify_idents import UnifiedRow
from unify_idents.engine_parsers.base_parser import __BaseParser
from pathlib import Path
import re
import csv
from decimal import Decimal, getcontext, ROUND_UP
import itertools
from loguru import logger

"""
1. Raw data location
2. Mass and m/z calculation
3. Remapping proteins
4. 

"""


class MSFragger3Parser(__BaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        if params is None:
            params = {}
        self.params = params
        self.input_file = input_file

        try:
            self.reader = csv.DictReader(open(input_file), delimiter="\t")
        except:
            self.reader = iter([])

        self.style = "msfragger_style_3"
        self.column_mapping = self.get_column_names()

        self.cols_to_add = [
            "Raw file location",
            "Spectrum Title",
            "uCalc m/z",
            "uCalc Mass",
            "Retention Time (s)",
            "Accuracy (ppm)",
            "Mass Difference",
            "Protein ID",
            "Sequence Start",
            "Sequence Stop",
            "Sequence Pre AA",
            "Sequence Post AA",
            "Enzyme Specificity",
            "Complies search criteria",
            "Conflicting uparam",
            "Search Engine",
        ]

        self.DICT_15N_DIFF = {
            "A": 0.997035,
            "C": 0.997035,
            "D": 0.997035,
            "E": 0.997035,
            "F": 0.997035,
            "G": 0.997035,
            "H": 2.991105,
            "I": 0.997035,
            "K": 1.99407,
            "L": 0.997035,
            "M": 0.997035,
            "N": 1.99407,
            "P": 0.997035,
            "Q": 1.99407,
            "R": 3.98814,
            "S": 0.997035,
            "T": 0.997035,
            "V": 0.997035,
            "W": 1.99407,
            "Y": 0.997035,
        }

        self.mass2mod = {str(d["mass"]): d["name"] for d in self.params["mods"]["opt"]}
        self.mass2mod.update(
            {str(d["mass"]): d["name"] for d in self.params["mods"]["fix"]}
        )
        self.mod_pattern = re.compile(r""":(?P<pos>[0-9]*$)""")
        self.mod_pattern_msfragger_term = re.compile(
            r""".-term\((?P<mass>[0-9]*\.[0-9]*)\)"""
        )
        self.mod_pattern_msfragger = re.compile(
            r"""(?P<pos>[0-9]*)(?P<aa>[A-Z])\((?P<mass>[0-9]*\.[0-9]*)\)"""
        )
        no_decimals = 4
        self.mass_format_string = "{{0:3.{0}f}}".format(no_decimals)
        self.mass_to_mod_combo = self.prepare_mass_to_mod()

    def __next__(self):
        line = next(self.reader)
        line = self._unify_row(line)
        return line

    @classmethod
    def file_matches_parser(cls, file):
        column_names = [
            "scannum",
            "peptide",
            "charge",
            "peptide_prev_aa",
            "peptide_next_aa",
            "protein",
            "modification_info",
            "retention_time",
            "precursor_neutral_mass",
            "calc_neutral_pep_mass",
            "hit_rank",
            "massdiff",
            "num_matched_ions",
            "tot_num_ions",
            "hyperscore",
            "nextscore",
            "num_tol_term",
            "num_missed_cleavages",
            "expectscore",
            "best_locs",
            "score_without_delta_mass",
            "best_score_with_delta_mass",
            "second_best_score_with_delta_mass",
            "delta_score",
        ]
        with open(file) as fin:
            headers = fin.readline()
            if set(headers.split()) == set(column_names):
                ret_val = True
            else:
                ret_val = False
        return ret_val

    def _unify_row(self, row):
        new_row = {}
        for unify_name, omssa_name in self.column_mapping.items():
            new_row[unify_name] = row[omssa_name]
        for col in self.cols_to_remove:
            del new_row[col]
        for col in self.cols_to_add:
            if col not in new_row:
                new_row[col] = ""
        new_row["Search Engine"] = "msfragger_3_0"
        new_row["Spectrum Title"] = "{file}.{specid}.{specid}.{charge}".format(
            file=self.params["Raw file location"].split(".")[0],
            specid=new_row["Spectrum ID"],
            charge=new_row["Charge"],
        )
        new_row["Raw file location"] = self.params["Raw file location"]
        new_row["Exp m/z"] = self.calc_mz(
            new_row["MSFragger:Precursor neutral mass (Da)"], new_row["Charge"]
        )
        new_row["Calc m/z"] = self.calc_mz(
            new_row["MSFragger:Neutral mass of peptide"], new_row["Charge"]
        )

        # TODO
        # think of all the labile mode, glycan and 15N stuff ...
        modstring = self.format_mods(new_row)

        new_row["Modifications"] = modstring

        new_row = self.general_fixes(new_row)
        return UnifiedRow(**new_row)

    def get_column_names(self):
        headers = self.param_mapper.get_default_params(style=self.style)[
            "header_translations"
        ]["translated_value"]
        return headers

    def prepare_mass_to_mod(self):

        self.fixed_mods = {}
        self.opt_mods = {}
        self.mod_dict = {}
        self.n_term_replacement = {
            "Ammonia-loss": None,
            "Trimethyl": None,
            "Gly->Val": None,
        }

        # 1.1 Create mod_dict with name as key and mass, aa, and pos as subkeys
        for mod_type in ["fix", "opt"]:
            for modification in self.params["mods"][mod_type]:
                aa = modification["aa"]
                pos = modification["pos"]
                name = modification["name"]
                if name not in self.mod_dict.keys():
                    self.mod_dict[name] = {
                        "mass": modification["mass"],
                        "aa": set(),
                        "pos": set(),
                    }
                self.mod_dict[name]["aa"].add(aa)

                # self.mod_dict[name]["aa"].add(pos)
                self.mod_dict[name]['pos'].add(pos)

                if "N-term" in pos:
                    self.n_term_replacement[name] = aa
                if mod_type == "fix":
                    self.fixed_mods[aa] = name
                    if aa == "C" and name == "Carbamidomethyl":
                        cam = True
                        self.mod_dict["Carbamidomethyl"]["aa"].add("U")
                        self.fixed_mods["U"] = "Carbamidomethyl"
                if mod_type == "opt":
                    self.opt_mods[aa] = name

        # rounding using Decimal (stdlib)
        getcontext().prec = 8
        getcontext().rounding = ROUND_UP
        # mod_dict_list = params['mods']['opt'] + params['mods']['fix']

        # 1.2 Add Mod 15N combinations to self.mod_dict
        # TODO get rid of use15N param
        use15N = self.params.get("15N", False)
        if use15N:
            aminoacids_2_check = set()
            for modname in self.mod_dict.keys():
                aminoacids_2_check |= self.mod_dict[modname]["aa"]
            additional_15N_modifications = []
            for aminoacid, N15_Diff in self.DICT_15N_DIFF.items():
                if aminoacid not in aminoacids_2_check:
                    continue
                if "_15N_{0}".format(aminoacid) in self.mod_dict.keys():
                    print(
                        """
                        Error in unify_csv
                        New mod_name already present in mod_dict'
                        This should never happen"""
                    )
                    sys.exit(1)
                self.mod_dict["_15N_{0}".format(aminoacid)] = {
                    "mass": N15_Diff,
                    "aa": set([aminoacid]),
                    "pos": set(["any"]),
                }

        # 2. Create all combinations of all possible mods, calc masses for these "merged mods" and save a mass to mod combo lookup
        mod_names = []
        for mod in sorted(list(self.mod_dict.keys())):
            mod_names.extend(itertools.repeat(mod, len(self.mod_dict[mod]["aa"])))
        mass_to_mod_combo = {}
        # we cover all combocs of mods
        for iter_length in range(2, len(mod_names) + 1):
            for name_combo in itertools.combinations(mod_names, iter_length):
                mass = 0
                for name in name_combo:
                    mass += Decimal(self.mod_dict[name]["mass"])
                rounded_mass = self.mass_format_string.format(mass)
                if rounded_mass not in mass_to_mod_combo.keys():
                    mass_to_mod_combo[rounded_mass] = set()
                mass_to_mod_combo[rounded_mass].add(name_combo)
        return mass_to_mod_combo

    def format_mods(self, row):
        final_mods = []
        for single_mod in row["Modifications"].split(", "):
            if single_mod.strip() == '':
                continue
            else:
                match = self.mod_pattern_msfragger.search(single_mod)
                pos = int(match.group("pos"))
                mass = float(match.group("mass"))
                aa = match.group("aa")
                mass_string = msfragger_mass = self.mass_format_string.format(
                    Decimal(mass)
                )

            if mass_string in self.mass_to_mod_combo.keys():
                explainable_combos = []
                for combo in self.mass_to_mod_combo[mass_string]:
                    for new_name in combo:
                        meta_mod_info = self.mod_dict[new_name]
                        real_pos = pos
                        if 'Prot-N-term' in self.mod_dict[new_name]['pos'] or 'N-term' in self.mod_dict[new_name]['pos']:
                            real_pos = 0 # should be 1 and is set to 0
                        final_mods.append(f'{new_name}:{real_pos}')
            else:
                # not a merged mod
                names = self.mod_mapper.appMass2name_list(mass, decimal_places=6)
                for new_name in names:
                    if new_name in self.mod_dict:
                        real_pos = pos
                        if 'Prot-N-term' in self.mod_dict[new_name]['pos'] or 'N-term' in self.mod_dict[new_name]['pos']:
                            real_pos = 0 # should be 1 and is set to 0
                        final_mods.append(f'{new_name}:{pos}')
                        break
        # breakpoint()
        # # 1.1 Extract all masses of mods and add to mod_list
        # # msfragger specific
        # ms_fragger_reformatted_mods = []
        # if row["Modifications"] == "M":
        #     # M stand for Modifications here, not Methionine
        #     row["Modifications"] = ""
        # else:
        #     mod_string = row["Modifications"]
        #     if mod_string == "":
        #         mod_list = []
        #     # elif "|" in mod_string:
        #     #     mod_list = []
        #     #     for single_mod in mod_string.split("|"):
        #     #         if single_mod in ["M", ""]:
        #     #             continue
        #     #         msfragger_pos, raw_msfragger_mass = single_mod.split("$")
        #     #         msfragger_mass = self.mass_format_string.format(
        #     #             # mass rounded as defined above
        #     #             Decimal(raw_msfragger_mass)
        #     #         )
        #     #         msfragger_pos = int(msfragger_pos)
        #     # mod_list.append((msfragger_mass, raw_msfragger_mass, msfragger_pos))
        #     else:
        #         mod_list = []
        #         for single_mod in mod_string.split(", "):
        #             if single_mod.startswith("N-term"):
        #                 msfragger_pos = 0
        #                 match = self.mod_pattern_msfragger_term.search(single_mod)
        #             elif single_mod.startswith("C-term"):
        #                 msfragger_pos = len(row["Sequence"])
        #             else:
        #                 match = self.mod_pattern_msfragger.search(single_mod)
        #                 msfragger_pos = int(match.group("pos"))
        #                 if msfragger_pos != 0:
        #                     msfragger_pos -= 1
        #             raw_msfragger_mass = float(match.group("mass"))
        #             msfragger_mass = self.mass_format_string.format(
        #                 # mass rounded as defined above
        #                 Decimal(raw_msfragger_mass)
        #             )
        #             mod_list.append((msfragger_mass, raw_msfragger_mass, msfragger_pos))
        #     # 1.2 For every mod mass, find out if its a combined mass and if yes, add all mods contributing to merged mod
        #     # pretty msfragger specific
        # ^^^^^^^^^^^^^^^^^^^^^^^ Done
        #     for single_mod in mod_list:
        #         (
        #             msfragger_mass,
        #             raw_msfragger_mass,
        #             msfragger_pos,
        #         ) = single_mod
        #         if msfragger_mass in self.mass_to_mod_combo.keys():
        #             explainable_combos = []
        #             for combo in self.mass_to_mod_combo[msfragger_mass]:
        #                 combo_explainable = set([True])
        #                 tmp_mods = []
        #                 for new_name in combo:
        #                     meta_mod_info = self.mod_dict[new_name]
        #                     single_mod_check = set([True])
        #                     # check aa
        #                     if (
        #                         "*" not in meta_mod_info["aa"]
        #                         and row["Sequence"][msfragger_pos]
        #                         not in meta_mod_info["aa"]
        #                     ):
        #                         single_mod_check.add(False)
        #                     # check pos
        #                     if "any" not in meta_mod_info["pos"]:
        #                         pos_to_check = set()
        #                         if (
        #                             "Prot-N-term" in meta_mod_info["pos"]
        #                             or "N-term" in meta_mod_info["pos"]
        #                         ):
        #                             pos_to_check.add(0)
        #                         elif (
        #                             "Prot-C-term" in meta_mod_info["pos"]
        #                             or "C-term" in meta_mod_info["pos"]
        #                         ):
        #                             pos_to_check.add(int(len(row["Sequence"])) - 1)
        #                         else:
        #                             pass
        #                         if pos_to_check != set():
        #                             if msfragger_pos not in pos_to_check:
        #                                 single_mod_check.add(False)

        #                     if all(single_mod_check):
        #                         pos_in_peptide_for_format_str = msfragger_pos + 1
        #                         # we keep mass here so that the
        #                         # correct name is added later in already
        #                         # existing code
        #                         tmp_mods.append(
        #                             "{0}:{1}".format(
        #                                 meta_mod_info["mass"],
        #                                 pos_in_peptide_for_format_str,
        #                             )
        #                         )
        #                     else:
        #                         combo_explainable.add(False)
        #                 if all(combo_explainable):
        #                     explainable_combos.append(tmp_mods)
        #             if len(explainable_combos) > 1:
        #                 print(
        #                     """
        #                     [ WARNING ] Multiple modification combinations possible
        #                     [ WARNING ] to explain reported modification mass
        #                     [ WARNING ] The following combination was chosen to continue:
        #                     [ WARNING ] {0}
        #                     """.format(
        #                         sorted(explainable_combos)[0],
        #                     )
        #                 )
        #             elif len(explainable_combos) == 1:
        #                 ms_fragger_reformatted_mods += sorted(explainable_combos)[0]
        #             else:
        #                 # no combos explainable
        #                 ms_fragger_reformatted_mods.append(
        #                     "{0}:{1}".format(raw_msfragger_mass, msfragger_pos + 1)
        #                 )
        #         else:
        #             # Not a merged mod
        #             ms_fragger_reformatted_mods.append(
        #                 "{0}:{1}".format(raw_msfragger_mass, msfragger_pos + 1)
        #             )
        # logger.debug(f"ms_fragger_reformatted_mods = {ms_fragger_reformatted_mods}")
        # # 2. Map mass to mod (move to base parser)
        # formatted_mods = []
        # for mod in ms_fragger_reformatted_mods:
        #     mass, pos = mod.split(":")
        #     mass, pos = float(mass), int(pos)
        #     names = self.mod_mapper.appMass2name_list(mass, decimal_places=6)
        #     logger.debug(names)
        #     for n in names:
        #         if n in self.mod_dict:
        #             if (
        #                 row["Sequence"][pos - 1] in self.mod_dict[n]["aa"]
        #                 or "*" in self.mod_dict[n]["aa"]
        #             ):
        #                 formatted_mods.append(f"{n}:{pos}")
        # row["Modifications"] = ";".join(formatted_mods)

        # 3. Check if pos is correct (including N-term/C-term)
        # for modification in row["Modifications"].split(";"):
        #     # could be moved to base parser
        #     ############################################################
        #     raw_modification = modification
        #     Nterm = False
        #     Cterm = False
        #     skip_mod = False
        #     if modification == "" or modification == "null":
        #         continue

        #     pos, mod = None, None
        #     match = self.mod_pattern.search(modification)
        #     pos = int(match.group("pos"))
        #     mod = modification[: match.start()]

        #     if pos <= 1:
        #         Nterm = True
        #         new_pos = 1
        #     elif pos > len(row["Sequence"]):
        #         Cterm = True
        #         new_pos = len(row["Sequence"])
        #     else:
        #         new_pos = pos
        #     aa = row["Sequence"][new_pos - 1].upper()
        #     ############################################################

        #     # TODO which of these cases actually occur in msfragger?
        #     # Could go to base_parser
        #     ############################################################
        #     if mod in self.mod_dict.keys():
        #         correct_mod = False
        #         if aa in self.mod_dict[mod]["aa"]:
        #             # everything is ok
        #             correct_mod = True
        #         elif Nterm or Cterm:
        #             if "*" in self.mod_dict[mod]["aa"]:
        #                 correct_mod = True
        #                 # still is ok
        #         assert (
        #             correct_mod is True
        #         ), """
        #                 A modification was reported for an aminoacid for which it was not defined
        #                 unify_csv cannot deal with this, please check your parameters and engine output
        #                 reported modification: {0} on {1}
        #                 modifications in parameters: {2}
        #                 """.format(
        #             mod, aa, params["mods"]
        #         )
        #     elif "unknown modification" == mod:
        #         modification_known = False
        #         if aa in opt_mods.keys():
        #             # fixed mods are corrected/added already
        #             modification = "{0}:{1}".format(opt_mods[aa], new_pos)
        #             modification_known = True
        #         assert (
        #             modification_known == True
        #         ), """
        #                 unify csv does not work for the given unknown modification for
        #                 {0} {1} aa: {2}
        #                 maybe an unknown modification with terminal position was given?
        #                 """.format(
        #             line_dict["Sequence"], modification, aa
        #         )
        #     else:
        #         float_mod = float(mod)
        #         masses_2_test = [float_mod]
        #         if use15N:
        #             substract_15N_diff = False
        #             if (
        #                 aa in fixed_mods.keys()
        #                 and "msgfplus" in search_engine
        #                 and pos != 0
        #             ):
        #                 substract_15N_diff = True
        #             if "msfragger" in search_engine and float_mod > 4:
        #                 # maximum 15N labeling is 3.988 Da (R)
        #                 substract_15N_diff = True
        #             if substract_15N_diff:
        #                 masses_2_test.append(float_mod - ursgal.ukb.DICT_15N_DIFF[aa])
        #         name_list = []
        #         for mass_2_test in masses_2_test:
        #             mass_buffer_key = mass_format_string.format(mass_2_test)
        #             # buffer increases speed massively...
        #             if mass_buffer_key not in app_mass_to_name_list_buffer.keys():
        #                 app_mass_to_name_list_buffer[
        #                     mass_buffer_key
        #                 ] = ursgal.GlobalUnimodMapper.appMass2name_list(
        #                     float(mass_buffer_key), decimal_places=no_decimals
        #                 )
        #             name_list += app_mass_to_name_list_buffer[mass_buffer_key]
        #         mapped_mod = False
        #         for name in name_list:
        #             if name in self.mod_dict.keys():
        #                 if aa in self.mod_dict[name]["aa"]:
        #                     modification = "{0}:{1}".format(name, new_pos)
        #                     mapped_mod = True
        #                 elif Nterm and "*" in self.mod_dict[name]["aa"]:
        #                     modification = "{0}:{1}".format(name, 0)
        #                     mapped_mod = True
        #                 else:
        #                     continue
        #             elif use15N and name in [
        #                 "Label:15N(1)",
        #                 "Label:15N(2)",
        #                 "Label:15N(3)",
        #                 "Label:15N(4)",
        #             ]:
        #                 mapped_mod = True
        #                 skip_mod = True
        #                 break
        #     ### Check if Mod is N-term, if yes, correct position (from 1 to 0)
        #     # could be moved to base parser
        #     if modification in tmp_mods:
        #         add_n_term_mod = False
        #         _mod, _pos = modification.split(":")
        #         _aa = line_dict["Sequence"][int(_pos) - 1]
        #         if _aa in self.mod_dict[_mod]["aa"]:
        #             if "N-term" in self.mod_dict[_mod]["aa"]:
        #                 add_n_term_mod = True
        #             elif "Prot-N-term" in self.mod_dict[_mod]["aa"]:
        #                 add_n_term_mod = True
        #         if add_n_term_mod:
        #             modification = modification.replace(
        #                 "{0}:1".format(_mod), "{0}:0".format(_mod)
        #             )
        #     if skip_mod is True:
        #         continue
        #     tmp_mods.append(modification)

        # # breakpoint()
        # # # if "msfragger" in search_engine:
        # # #     org_mass_diff = row["Mass Difference"]
        # # #     tmp_mass_diff.append("{0}:n".format(org_mass_diff))

        # # breakpoint()
        # formatted_mods = ";".join(tmp_mods)
        return ";".join(final_mods)
