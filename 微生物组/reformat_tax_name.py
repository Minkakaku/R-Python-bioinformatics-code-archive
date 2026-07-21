# Python script to reformat taxonomy name output by QIIME2 in a more human-readable format
# By Yewei on 2024.10.28

# Example:
# "d__Eukaryota;p__Ascomycota;__" --> "d__Eukaryota;p__Ascomycota;c__unclassified_Ascomycota"
# "Unassigned;__;__" --> "Unassigned;Unassigned;Unassigned"


def reformat_tax_name(raw, mode, mode2=""):
    import pandas as pd
    tax_list = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    if mode == "16S" or mode == "18S":
        tax_short_list = ["d", "p", "c", "o", "f", "g", "s"]
    elif mode == "ITS":
        tax_short_list = ["k", "p", "c", "o", "f", "g", "s"]
    # If taxonomy name in one level starts with these terms, it will be converted to "unclassfied_{upperLevel}"
    target_list = ["uncultured", "unidentified", "gut_metagenome", "metagenome", "mouse_gut", "bacterium_enrichment",
                   "human_gut", "Incertae_Sedis", "Unknown"]
    if mode2 == "ASV":
        level_count = 7
    else:
        level_count = len(raw.split(';'))
    # Raw taxonomy name will be converted to two parts: reformatted version ("Taxonomy") and specific level ("{Level}")
    raw_list = raw.split(";")
    if mode2 == "ASV" and len(raw_list) < 7:
        raw_list = raw_list + ["__"] * (7 - len(raw_list))
    new_list = []
    if raw_list[0] == "Unassigned":
        tax_new = ';'.join(["Unassigned"] * level_count)
    else:
        for i in range(level_count):
            raw_item = raw_list[i].split("__")
            if (i > 0) and (raw_item[0] == "" or raw_item[1].startswith(tuple(target_list))):
                if new_list[i - 1].split("__")[1].startswith("unclassified"):
                    # If last level has already been reformatted (starting with "unclassified"), just keep it
                    new = "{}__{}".format(tax_short_list[i], new_list[i - 1].split("__")[1])
                else:
                    new = "{}__unclassified_{}".format(tax_short_list[i], new_list[i - 1].split("__")[1])
            elif (i > 0) and (raw_item[1] == new_list[i - 1].split("__")[1]):
                # If this level has the same name as last level, add "unclassified_" as prefix
                new = "{}__unclassified_{}".format(tax_short_list[i], new_list[i - 1].split("__")[1])
            elif (i > 0) and (raw_item[1] == new_list[i - 1].split("__")[1].replace("unclassified_", "")):
                # If this level has the same name as last level, add "unclassified_" as prefix
                new = "{}__{}".format(tax_short_list[i], new_list[i - 1].split("__")[1])
            else:
                new = raw_list[i]
            new_list.append(new)
        tax_new = ';'.join(new_list)
    return tax_new

