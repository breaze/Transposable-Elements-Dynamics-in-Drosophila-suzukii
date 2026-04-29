#!/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


classification_dict = {"CLASSII/TIR/ACADEM-1": ["Academ", "DNA/Academ-1", "CLASSII/TIR/ACADEM-1", "DNA/Academ-2", "DNA/Academ-H"],
					   "CLASSII/TIR/TIR": ["DNA", "CLASSII/TIR/TIR", "CLASSII/TIR", "DNA/DNA"],
                       "CLASSII/TIR/CMC": ["CLASSII/TIR/CMC", "DNA/CMC-Chapaev", "DNA/CMC-Chapaev-3", "DNA/CMC-EnSpm", "DNA/CMC-Mirage", "DNA/CMC-Transib", "DNA/CMC"],
                       "CLASSII/TIR/TRANSIB": ["CLASSII/TIR/TRANSIB", "Transib", "DNA/Transib"],
					   "CLASSII/CRYPTON/CRYPTON": ["CryptonA", "Crypton", "CLASSII/CRYPTON/CRYPTON", "CLASSII/CRYPTON", "DNA/Crypton-A", "DNA/Crypton-H", "DNA/Crypton-V"],
					   "CLASSII/TIR/DADA": ["CLASSII/TIR/DADA", "DNA/Dada", "Dada"],
					   "CLASSII/TIR/GINGER": ["CLASSII/TIR/GINGER", "DNA/Ginger-1", "DNA/Ginger-2"],
					   "CLASSII/TIR/IS3EU": ["CLASSII/TIR/IS3EU", "DNA/IS3EU", "IS3EU"],
					   "CLASSII/TIR/KOLOBOK": ["Kolobok", "CLASSII/TIR/KOLOBOK", "DNA/Kolobok", "DNA/Kolobok-E", "DNA/Kolobok-Hydra", "DNA/Kolobok-T2", "DNA/Kolobok-H"],
					   "CLASSII/TIR/MULE": ["MuDR", "CLASSII/TIR/MULE", "DNA/MULE", "DNA/MULE-MuDR", "DNA/MULE-NOF", "TIR/MuDR", "DNA/MuDR", "CLASSII/TIR/MUDR", "CLASSII/TIR/MUTATOR", "DNA/MULE-F"],
					   "CLASSII/MAVERICK/MAVERICK": ["Polinton", "CLASSII/MAVERICK/MAVERICK", "DNA/Maverick", "CLASSII/MAVERICK"],
					   "CLASSII/TIR/MERLIN": ["Merlin", "CLASSII/TIR/MERLIN", "DNA/Merlin"],
					   "CLASSII/TIR/P": ["P", "CLASSII/TIR/P", "DNA/P", "TIR/P"],
					   "CLASSII/TIR/P-FUNGI": ["DNA/P-Fungi", "CLASSII/TIR/P-FUNGI"],
					   "CLASSII/TIR/CACTA": ["CLASSII/TIR/CACTA", "TIR/EnSpm/CACTA", "EnSpm", "DNA/EnSpm/CACTA", "EnSpm/CACTA"],
					   "CLASSII/TIR/PIF": ["ISL2EU", "CLASSII/TIR/PIF", "DNA/PIF", "DNA/PIF-ISL2EU", "DNA/PIF-Spy", "DNA/PIF-HarbS"],
					   "CLASSII/TIR/PIFHARBINGER": ["Harbinger", "CLASSII/TIR/PIFHARBINGER", "DNA/PIF-Harbinger", "TIR/Harbinger", "DNA/Harbinger"],
					   "CLASSII/TIR/PIGGYBAC": ["piggyBac", "CLASSII/TIR/PIGGYBAC", x"DNA/PiggyBac", "DNA/PiggyBac-X"],
					   "CLASSII/TIR/SOLA": ["Sola3", "Sola1", "Sola2", "CLASSII/TIR/SOLA", "DNA/Sola-1", "DNA/Sola-2", "DNA/Sola-3"],
					   "CLASSII/TIR/TC1MARINER": ["Mariner/Tc1", "CLASSII/TIR/TC1MARINER", "DNA/TcMar", "DNA/TcMar-Fot1", "DNA/TcMar-ISRm11", "DNA/TcMar-Mariner", "DNA/TcMar-Pogo", "DNA/TcMar-Tc1", "DNA/TcMar-Tc2", "DNA/TcMar-Tc4", "DNA/TcMar-Tigger", "DNA/TcMar-m44", "TIR/Mariner/Tc1", "CLASSII/TIR/TC1-MARINER", "DNA/TcMar-Ant1", "DNA/TcMar-Stowaway", "DNA/TcMar-Cweed", "DNA/TcMar-Sagan", "DNA/TcMar-Sagan", "DNA/Tc1-Mariner"],
					   "CLASSII/TIR/ZATOR": ["CLASSII/TIR/ZATOR", "DNA/Zator"],
					   "CLASSII/TIR/ZISUPTON": ["Zisupton", "CLASSII/TIR/ZISUPTON", "DNA/Zisupton"],
					   "CLASSII/TIR/HAT": ["hAT", "CLASSII/TIR/HAT", "DNA/hAT", "DNA/hAT-Ac", "DNA/hAT-Blackjack", "DNA/hAT-Charlie", "DNA/hAT-Tag1", "DNA/hAT-Tip100", "DNA/hAT-hAT19", "DNA/hAT-hAT5", "DNA/hAT-hAT6", "DNA/hAT-hATm", "DNA/hAT-hATw", "DNA/hAT-hATx", "DNA/hAT-hobo", "TIR/hAT", "DNA/hAT-Pegasus", "DNA/hAT-hAT1", "DNA/hAT-Restless"],
					   "CLASSII/HELITRON/HELITRON": ["Helitron", "CLASSII/HELITRON/HELITRON", "CLASSII/HELITRON", "Helitron/Helitron", "RC/Helitron", "RC/Helentron"],
					   "CLASSI/LINE/LINE": ["Non-LTR", "CLASSI/LINE/LINE", "LINE", "CLASSI/LINE", "LINE/LINE"],
					   "CLASSI/LINE/CR1": ["CLASSI/LINE/CR1", "LINE/CR1", "LINE/CR1-Zenon", "CR1"],
					   "CLASSI/LINE/DONG-R4": ["CLASSI/LINE/DONG-R4", "LINE/Dong-R4"],
					   "CLASSI/LINE/I": ["CLASSI/LINE/I", "LINE/I", "Nimb", "I"],
					   "CLASSI/LINE/JOCKEY": ["CLASSI/LINE/JOCKEY", "LINE/I-Jockey", "LINE/Jockey"],
					   "CLASSI/LINE/L1": ["CLASSI/LINE/L1", "LINE/L1", "LINE/L1-Tx1", "L1", "Tx1", "LINE/L1-Zorro"],
					   "CLASSI/LINE/L2": ["CLASSI/LINE/L2", "LINE/L2", "L2"],
					   "CLASSI/LINE/DAPHNE": ["CLASSI/LINE/DAPHNE", "Daphne"],
					   "CLASSI/LINE/VINGI": ["CLASSI/LINE/VINGI", "Vingi"],
					   "CLASSI/LINE/CRACK": ["CLASSI/LINE/CRACK", "Crack"],
					   "CLASSI/LINE/HERO": ["CLASSI/LINE/HERO", "Hero"],
					   "CLASSI/LINE/NESL": ["CLASSI/LINE/NESL", "NeSL"],
					   "CLASSI/LINE/PROTO1": ["CLASSI/LINE/PROTO1", "LINE/Proto1", "Proto1"],
					   "CLASSI/LINE/PROTO2": ["CLASSI/LINE/PROTO2", "LINE/Proto2", "Proto2"],
					   "CLASSI/LINE/R1": ["CLASSI/LINE/R1", "LINE/R1", "LINE/R1-LOA", "R1"],
					   "CLASSI/LINE/R2": ["CLASSI/LINE/R2", "LINE/R2", "LINE/R2-Hero", "LINE/R2-NeSL", "R2"],
					   "CLASSI/LINE/RTE": ["RTE", "CLASSI/LINE/RTE", "LINE/RTE", "LINE/RTE-BovB", "LINE/RTE-RTE", "LINE/RTE-X", "LINE/RTE-ORTE"],
					   "CLASSI/LINE/RTEX": ["CLASSI/LINE/RTEX", "RTEX"],
					   "CLASSI/LINE/REX": ["CLASSI/LINE/REX", "Rex1", "Rex2", "Rex3", "Rex6"],
					   "CLASSI/LINE/CRE": ["CLASSI/LINE/CRE", "CRE", "LINE/CRE-Odin", "LINE/CRE"],
					   "CLASSI/LINE/RANDI": ["CLASSI/LINE/RANDI", "RandI"],
					   "CLASSI/LINE/TAD1": ["CLASSI/LINE/Tad1", "LINE/Tad1"],
					   "CLASSI/LINE/REX-BABAR": ["CLASSI/LINE/REX-BABAR", "LINE/Rex-Babar"],
					   "CLASSI/LINE/DECEIVER": ["CLASSI/LINE/DECEIVER", "LINE/Deceiver"],
					   "CLASSI/LINE/DUALEN": ["CLASSI/LINE/DUALEN", "LINE/Dualen"],
					   "CLASSI/LINE/LOA": ["CLASSI/LINE/LOA", "LINE/LOA"],
					   "CLASSI/LTR/LTR": ["CLASSI/LTR/LTR", "LTR", "CLASSI/LTR", "LTR/LTR"],
					   "CLASSI/LTR/CAULIMOVIRUS": ["CLASSI/LTR/CAULIMOVIRUS", "Caulimovirus", "LTR/Caulimovirus", "CLASSI/LTR/Caulimovirus"],
					   "CLASSI/LTR/COPIA": ["CLASSI/LTR/COPIA", "LTR/Copia", "Copia"],
					   "CLASSI/DIRS/DIRS": ["DIRS", "CLASSI/DIRS/DIRS", "LTR/DIRS"],
					   "CLASSI/DIRS/NGARO": ["CLASSI/DIRS/NGARO", "LTR/Ngaro", "CLASSI/DIRS/Ngaro"],
					   "CLASSI/DIRS/Unknown": ["CLASSI/DIRS/Unknown"],
					   "CLASSI/LTR/ERV": ["ERV1", "ERV2", "ERV3", "Endogenous", "Endogenous Retrovirus", "CLASSI/LTR/ERV", "LTR/ERV", "LTR/ERV-Foamy", "LTR/ERV1", "LTR/ERV4", "LTR/ERV-Foamy", "LTR/ERV-MaLR", "LTR/ERVK", "LTR/ERVL", "LTR/ERVL-MaLR", "LTR/Endogenous"],
					   "CLASSI/LTR/GYPSY": ["CLASSI/LTR/GYPSY", "LTR/Gypsy", "Gypsy"],
					   "CLASSI/LTR/BELPAO": ["CLASSI/LTR/BELPAO", "LTR/Pao", "BEL", "CLASSI/LTR/BEL-PAO", "LTR/Bel-Pao"],
					   "CLASSI/PLE/PLE": ["Penelope/Poseidon", "Neptune", "Coprina", "Naiad/Chlamys", "CLASSI/PLE/PLE", "Penelope", "PLE", "PLE/Chlamys", "PLE/Hydra", "PLE/Naiad", "CLASSI/PLE"],
					   "CLASSI/RETROPOSON/RETROPOSON": ["CLASSI/RETROPOSON/RETROPOSON", "Retroposon", "Retroposon/L1-dep", "Retroposon/L1-derived", "Retroposon/RTE-derived", "Retroposon/SVA", "Retroposon/sno", "Retroposon/sno-L1"],
					   "CLASSI/SINE/SINE": ["CLASSI/SINE", "CLASSI/SINE/SINE", "SINE", "SINE/SINE"],
					   "CLASSI/SINE/5S": ["CLASSI/SINE/5S", "SINE/5S", "SINE/5S-Deu-L2"],
					   "CLASSI/SINE/7S": ["CLASSI/SINE/7S", "SINE/7SL"],
					   "CLASSI/SINE/ALU": ["CLASSI/SINE/ALU", "SINE/Alu"],
					   "CLASSI/SINE/B2": ["CLASSI/SINE/B2", "SINE/B2"],
					   "CLASSI/SINE/B4": ["CLASSI/SINE/B4", "SINE/B4"],
					   "CLASSI/SINE/ERV1": ["CLASSI/SINE/ERV1", "SINE/ERV1"],
					   "CLASSI/SINE/ID": ["CLASSI/SINE/ID", "SINE/ID"],
					   "CLASSI/SINE/MIR": ["CLASSI/SINE/MIR", "SINE/MIR"],
					   "CLASSI/SINE/U": ["CLASSI/SINE/U", "SINE/U", "SINE/U-L1"],
					   "CLASSI/SINE/tRNA": ["SINE2/tRNA","CLASSI/SINE/tRNA", "SINE/tRNA", "SINE/tRNA-Core", "SINE/tRNA-Core-RTE", "SINE/tRNA-Deu", "SINE/tRNA-Deu-CR1", "SINE/tRNA-Deu-L2", "SINE/tRNA-Deu-RTE", "SINE/tRNA-L1", "SINE/tRNA-Meta", "SINE/tRNA-RTE", "SINE/tRNA-V", "SINE/tRNA-V-RTE", "SINE/SINE2/tRNA", "SINE/tRNA-CR1", "SINE/tRNA-V-CR1"],
					   "CLASSI/LTR/LARD": ["CLASSI/LTR/LARD"],
					   "CLASSI/LTR/TRIM": ["CLASSI/LTR/TRIM"],
					   "CLASSI/LTR/TR-GAG": ["CLASSI/LTR/TR-GAG"],
					   "CLASSII/MITE/MITE": ["miniature", "CLASSII/MITE/MITE", "CLASSII/MITE", "DNA/MITE"],
					   "CLASSII/REPLITON/REPLITON" : ["CLASSII/REPLITON/REPLITON", "Replitron"],
					   "CLASSII/TIR/NOVOSIB": ["CLASSII/TIR/NOVOSIB", "Novosib"],
					   "UNCLASSIFIED/UNCLASSIFIED/UNCLASSIFIED" : ["UNCLASSIFIED"],
					   "EXCLUDED": ["MSAT", "Transposable", "Unknown", "snRNA", "Satellite", "SAT", "Simple_Repeat", "Caulimoviridae", "Multicopy_gene", "Multicopy", "REP-10_OS", "RCH2", "Nonautonomous", "KRISPIE", "REP-1_OS", "ATREP18", "Rep-1_AT", "ATREP19", "tRNA"]
					   }

def standarization(lib_filename, oneCode):
	standard_seqs = []
	#species = os.path.basename(lib_filename).replace(".fasta", "").replace(".fa", "").replace("_", " ")
	for transposon in SeqIO.parse(lib_filename, "fasta"):
		classification = transposon.id.split("#")[1].split(" ")[0]
		species = " ".join(transposon.description.split(" ")[1:])
		found = 0
		if classification not in classification_dict["EXCLUDED"]:
			for new_class in classification_dict.keys():
				if classification in classification_dict[new_class]:
					found = 1
					transposon.id = transposon.id.split("#")[0] + "#" + new_class
					transposon.description = species
					standard_seqs.append(transposon)
					break

	
		if found == 0 and classification not in classification_dict["EXCLUDED"]:
			print("The TE named " + transposon.id + " from the species "+ species +" was not found !!!!")

	if oneCode:
		print("Doing the OneCode parse ...")
		standard_seqs = oneCode_parser(standard_seqs)
	return standard_seqs

def join_LTR_INT(lib_filename, standard_seqs):
	result_file = lib_filename+"_std.fa"
	species = os.path.basename(lib_filename).replace(".fasta", "").replace(".fa", "").replace("_", " ")
	dicc_LTR = {}
	dicc_int = {}

	for seq in standard_seqs:
		if "CLASSI/LTR" in seq.id.split("#")[1]:
			seqname = str(seq.id).split("#")[0]
			if "_INT" in seqname:
				real_seqname = seqname.replace("_INT", "")
				dicc_int[real_seqname.upper()] = str(seq.seq)
			elif "_int" in seqname:
				real_seqname = seqname.replace("_int", "")
				dicc_int[real_seqname.upper()] = str(seq.seq)
			elif "_LTR" in seqname:
				real_seqname = seqname.replace("_LTR", "")
				dicc_LTR[real_seqname.upper()] = str(seq.seq)

			# experimental case
			elif "-LTR" in seqname:
				real_seqname = seqname.replace("-LTR", "")
				dicc_LTR[real_seqname.upper()] = str(seq.seq)
			elif "_I" in seqname:
				real_seqname = seqname.replace("_I", "")
				dicc_int[real_seqname.upper()] = str(seq.seq)
			elif "-I" in seqname:
				real_seqname = seqname.replace("-I", "")
				dicc_int[real_seqname.upper()] = str(seq.seq)

	seq_list = []
	already_done = []
	for seq in standard_seqs:
		seqname = str(seq.id).split("#")[0].replace("_INT", "").replace("_int", "").replace("_LTR", "").replace("-LTR", "").replace("_I", "").replace("-I", "")
		classification = str(seq.id).split("#")[1]
		if seqname.upper() in dicc_LTR.keys() and seqname.upper() in dicc_int.keys():
			if seqname.upper() not in already_done:
				newseq = SeqRecord(Seq(dicc_LTR[seqname.upper()]+dicc_int[seqname.upper()]+dicc_LTR[seqname.upper()]), id=seqname+"#"+classification)
				newseq.description = species
				seq_list.append(newseq)
				already_done.append(seqname.upper())
		else:
			#seq.description = species
			seq_list.append(seq)
	SeqIO.write(seq_list, result_file, "fasta")

def oneCode_parser(input_seqs):
	seq_list = []
	oneCode_dict = {"HELITRON/HELITRON": "RC/Helitron", "MITE/MITE": "DNA/MITE", "MAVERICK/MAVERICK": "DNA/Maverick", "DIRS/NGARO": "LTR/NGARO", 
	"CRYPTON/CRYPTON": "DNA/Crypton", "PLE/PLE": "LINE/PLE", "DIRS/DIRS": "LTR/DIRS", "UNCLASSIFIED/UNCLASSIFIED": "Unknown"}

	for record in input_seqs:
		header_class = "/".join(record.id.split("#")[1].split(" ")[0].split("/")[1:]).replace("TIR", "DNA")
		header_speci = record.description

		if header_class in oneCode_dict.keys():
			header_class = oneCode_dict[header_class]

		record.id = record.id.split("#")[0] + "#" + header_class
		record.description = header_speci
		seq_list.append(record)

	return seq_list

lib_filename = sys.argv[1]
oneCode = sys.argv[2]
if oneCode.upper() in ["T", "TRUE"]:
	oneCode = True
else:
	print("INFO: OneCode output was set to false")
	oneCode = False
stand_seqs = standarization(lib_filename, oneCode)
join_LTR_INT(lib_filename, stand_seqs)
