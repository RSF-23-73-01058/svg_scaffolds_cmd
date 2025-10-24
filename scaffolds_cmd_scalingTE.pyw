import argparse
import base64
import csv
import io
import json
import re
import signal
import sys
import time
import zlib
from collections import Counter
from func_timeout import func_timeout, FunctionTimedOut
from operator import itemgetter
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import Geometry

##Input
##Input
data_input = sys.stdin.read()



def process_input(data):
	#Parse the input
	data_base = json.loads(data)
	##Settings
	#Adjustable
	non_generic = 1 if data_base['meta']['scaffold_type'] == "non-generic" else 0;
	max_range = int(data_base['meta']['max_range'])
	#Soldered
	structure_limit = 10000
	scaffold_limit = 500
	square_size = 200
	color = (1, 1, 1, 1)
	## Calculate scaffolds
	#	1.		Get the structures in SMILES from TAB and calculate scaffolds
	#	2.		Get mol
	#	3.1.	Count if mol is None
	#	3.		Get number of components
	#	4.		Get the largest component
	#	5.		Get the MW
	#	6.		Check MW (1800 for now)
	#	6.1.	Count if MW > 1800
	#	7.		Get number of CC bonds
	#	8.		Check number of CC bonds (1 for now)
	#	8.1.	Count if CC < 1
	#	9.		Remove stereo, Hs
	#	10.		Get scaff
	#	11.		Repeating SMILES: Count the number of SMILES repeating other SMILES: if there are 2 exact matches -> 1 repeat; 3 exact matches -> 2 repeats, 2 & 3 distinct exact matches -> 3 repeats
	##Some vars to store the intermediate results
	scaff_dict = {}
	scaffold_counted = []
	#scaff_dict['scaff_smiles'][0] -> number of structures;
	#scaff_dict['scaff_smiles'][1] -> list of structure_ids;
	smiles_list = []
	data_to_transform = []
	data_to_json = []
	acyclic_ids = []
	##Vars for statistics' accumulation
	failed_mol = 0
	failed_components = 0
	failed_mw = 0
	failed_bond = 0
	failed_stereo = 0
	failed_hs = 0
	failed_smiles = 0
	failed_scaff = 0
	failed_scaff_smiles = 0
	failed_scaff_mol = 0
	n_structures = 0
	interval_1 = 0;
	interval_2 = 0;
	interval_3 = 0;
	interval_4 = 0;
	interval_5 = 0;
	interval_6 = 0;
	interval_7 = 0;
	interval_8 = 0;
	interval_9 = 0;
	interval_10 = 0;
	interval_11 = 0;
	interval_12 = 0;
	interval_13 = 0;
	##Loop through the input structures
	for row in data_base['data']:
		try:
			current_mol = Chem.MolFromSmiles(base64.b64decode(row[1]))
		except:
			failed_mol += 1
			continue
		current_structure_id = row[0]
		n_structures += 1
		if n_structures > structure_limit:
			n_structures -= 1
			break
		try:
			n_components = len(Chem.rdmolops.GetMolFrags(current_mol))
		except:
			failed_components += 1
			continue
		if n_components > 1:
			current_mol = rdMolStandardize.LargestFragmentChooser().choose(current_mol)
		if Descriptors.MolWt(current_mol) > 1800:
			failed_mw += 1
			continue
		n_cc_bonds = 0
		n_bonds = 0
		for bond in current_mol.GetBonds():
			n_bonds += 1
			anum1 = bond.GetBeginAtom().GetAtomicNum()
			anum2 = bond.GetEndAtom().GetAtomicNum()
			if anum1 == 6 and anum2 == 6:
				n_cc_bonds += 1
		if n_cc_bonds < 1:
			failed_bond += 1
			continue
		n_cc_bonds = 0
		try:
			dummy_plain = Chem.RemoveStereochemistry(current_mol)
		except:
			failed_stereo += 1
			continue
		try:
			current_mol = Chem.RemoveHs(current_mol)
		except:
			failed_hs += 1
			continue
		try:	
			smiles = Chem.MolToSmiles(current_mol)
		except:	
			failed_smiles += 1
			continue
		try:
			scaff = MurckoScaffold.GetScaffoldForMol(current_mol)
		except:
			failed_scaff += 1
			continue
		if non_generic == 0:
			try:
				scaff = MurckoScaffold.MakeScaffoldGeneric(scaff)
			except:
				failed_scaff += 1
				continue
			try:
				scaff = MurckoScaffold.GetScaffoldForMol(scaff)
			except:
				failed_scaff += 1
				continue
		try:		
			scaff_smiles = Chem.MolToSmiles(scaff)
		except:
			failed_scaff_smiles += 1
			continue
		try:	
			check_scaff_byMol = Chem.MolFromSmiles(scaff_smiles)
		except:
			failed_scaff_mol += 1
			continue
		if scaff_smiles == '': scaff_smiles = '*'
		#Check if this scaffold has been seen already
		if scaff_smiles in scaff_dict:
			scaff_dict[scaff_smiles][0] = scaff_dict[scaff_smiles][0] + 1
			scaff_dict[scaff_smiles][1].append(current_structure_id)
		else:
			scaff_dict[scaff_smiles] = []
			scaff_dict[scaff_smiles].append(1)
			scaff_dict[scaff_smiles].append([])
			scaff_dict[scaff_smiles][1].append(current_structure_id)
		smiles_list.append(smiles)
	##Gather meta data
	#Number of processed structures
	n_structures = n_structures - failed_mol - failed_components - failed_mw - failed_bond - failed_stereo - failed_hs - failed_smiles - failed_scaff - failed_scaff_smiles - failed_scaff_mol
	if n_structures == 0:
		output_json = json.dumps({'status': 'EMPTY', 'meta': '', 'counts': '', 'depicted_description': '', 'data': "b'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'"})
		print(output_json)
	#Check for repeating SMILES
	smiles_counted = Counter(smiles_list).most_common()
	smiles_repeats = 0
	for smiles in smiles_counted:
		if smiles[1] < 2:
			break
		smiles_repeats += smiles[1] - 1
	meta_data = {'scaffold_type': data_base['meta']['scaffold_type'], 'successfully_processed': n_structures, 'repeats': smiles_repeats, 'failed_parse': failed_mol, 'failed_components': failed_components, 'failed_mw': failed_mw, 'failed_bond': failed_bond, 'failed_stereo': failed_stereo, 'failed_hs': failed_hs, 'failed_smiles': failed_smiles, 'failed_scaff': failed_scaff, 'failed_scaff_smiles': failed_scaff_smiles, 'failed_scaff_mol': failed_scaff_mol}
	#Count scaffolds
	#Sort the dict, SEE: https://stackoverflow.com/questions/1217251/sorting-a-dictionary-with-lists-as-values-according-to-an-element-from-the-list
	scaff_dict = dict(sorted(scaff_dict.items(), key=lambda e: e[1][0], reverse=True))
	for key, value in scaff_dict.items():
		temp = [key, value[0], value[1]]
		scaffold_counted.append(temp)
	n_scaff = len(scaffold_counted)
	acyclic_count = 0
	unknown_count = 0
	undepictable_count = 0
	#Calculate score
	max_p = scaffold_counted[0][1] / n_structures
	min_p = scaffold_counted[-1][1] / n_structures
	#Calculate score for the most abundant scaffolds
	if len(scaffold_counted) > scaffold_limit: 
		max_p_top = scaffold_counted[0][1] / n_structures
		min_p_top = scaffold_counted[scaffold_limit-1][1] / n_structures
	else:
		max_p_top = max_p
		min_p_top = min_p
	#Get the data for barchart
	for scaff_count in scaffold_counted:
		if scaff_count[0] == '*':
			acyclic_count = scaff_count[1] 
			continue
		if scaff_count[0] == '[*-]':
			unknown_count = scaff_count[1] 
			continue
		if scaff_count[0] == '[*+]':
			undepictable_count = scaff_count[1] 
			continue
		match scaff_count[1]:
			case _ if scaff_count[1] > 1000:
				interval_13 += 1
			case _ if scaff_count[1] > 100 and scaff_count[1] <= 1000:
				interval_12 += 1
			case _ if scaff_count[1] > 10 and scaff_count[1] <= 100:
				interval_11 += 1
			case _ if scaff_count[1] == 10:
				interval_10 += 1
			case _ if scaff_count[1] == 9:
				interval_9 += 1
			case _ if scaff_count[1] == 8:
				interval_8 += 1
			case _ if scaff_count[1] == 7:
				interval_7 += 1
			case _ if scaff_count[1] == 6:
				interval_6 += 1
			case _ if scaff_count[1] == 5:
				interval_5 += 1
			case _ if scaff_count[1] == 4:
				interval_4 += 1
			case _ if scaff_count[1] == 3:
				interval_3 += 1
			case _ if scaff_count[1] == 2:
				interval_2 += 1
			case _ if scaff_count[1] == 1:
				interval_1 += 1
	struct_per_scaff = n_structures / n_scaff
	counts = { 'total_amount': n_scaff, 'acyclic_count': acyclic_count, 'struct_per_scaff': struct_per_scaff, 'bar_counts': {'interval_1':interval_1/n_scaff, 'interval_2':interval_2/n_scaff, 'interval_3':interval_3/n_scaff,
	'interval_4':interval_4/n_scaff, 'interval_5':interval_5/n_scaff, 'interval_6':interval_6/n_scaff, 'interval_7':interval_7/n_scaff, 'interval_8':interval_8/n_scaff, 'interval_9':interval_9/n_scaff,
	'interval_10':interval_10/n_scaff, 'interval_11':interval_11/n_scaff, 'interval_12':interval_12/n_scaff, 'interval_13':interval_13/n_scaff} }
	#Determine the fraction of the selected scaffolds
	if scaffold_limit >= len(scaffold_counted):
		fraction_selected_scaffold = 1
	else:
		fraction_selected_scaffold = scaffold_limit / len(scaffold_counted)
	#Choose scaffolds to work with
	if fraction_selected_scaffold < 1:
		scaffold_counted = scaffold_counted[0:scaffold_limit]
		#Update the number of structures
		n_structures_selected = sum(list(zip(*scaffold_counted))[1])
	else:
		n_structures_selected = n_structures
	#What fraction of structures is described by the selected scaffolds?
	fraction_selected_structure = n_structures_selected / n_structures
	#Describe the resulting scaffolds' subset
	selected_description = {'scaffold_limit': scaffold_limit, 'scaffolds_selected': len(scaffold_counted), 'scaffold_fraction': fraction_selected_scaffold, 'structure_fraction_forSelected':  fraction_selected_structure}
	#Some vars to process the scaffolds
	scaffold_id = 0
	#Should the scores be mapped on the interval [1, max_range]?
	size_adjuster = 1 / min_p_top
	size_mapping = 0
	if max_p_top / min_p_top > max_range : size_mapping = 1
	##Do the actual drawing and scaling for each scaffold
	for scaff_count in scaffold_counted:
		scaffold_status = 1
		scaffold_id += 1
		scaffold_id_str = "scaff_" + str(scaffold_id) 
		if size_mapping == 1:
			score = ((scaff_count[1] / n_structures) - min_p_top) * ((max_range-1)/(max_p_top - min_p_top)) + 1
		else:
			score = (scaff_count[1] / n_structures) * size_adjuster
		scaff_mol = Chem.MolFromSmiles(scaff_count[0])
		#Process the case where scaffold is not present
		if scaff_count[0] == '*': scaff_mol.GetAtoms()[0].SetProp('atomLabel', "ACYCLIC")
		dummy_stereo = AllChem.Compute2DCoords(scaff_mol, useRingTemplates=True)
		#Enumerate atoms
		atoms = []
		for a in scaff_mol.GetAtoms():
			atoms.append(a.GetIdx())
		#Enumerate rings
		rings = ()
		rInfo = scaff_mol.GetRingInfo()
		rings = rInfo.AtomRings()
		#Prepare colors for atoms
		atom_cols = {}
		for i, at in enumerate(atoms):
			atom_cols[at] = color
		#Draw
		drawer = rdMolDraw2D.MolDraw2DSVG(-1,-1)
		drawer.drawOptions().padding = 0
		#Process things having structure
		if scaff_count[0] != '*':
			drawer.drawOptions().fillHighlights=True
			drawer.drawOptions().setHighlightColour((color))
			drawer.drawOptions().highlightBondWidthMultiplier = 4
			#Draw the molecule
			drawer.DrawMolecule(scaff_mol, highlightAtoms=atoms, highlightAtomColors=atom_cols)
			#Fill the rings
			if rings:
				drawer.ClearDrawing()
				conf = scaff_mol.GetConformer()
				for aring in rings:
					ps = []
					for ratom in aring:
						pos = Geometry.Point2D(conf.GetAtomPosition(ratom))
						ps.append(pos)
					#Actually fill them
					drawer.SetFillPolys(True)
					drawer.SetColour(color)
					drawer.DrawPolygon(ps)
			drawer.FinishDrawing()
		else:
			#Draw the molecule
			drawer.DrawMolecule(scaff_mol)
			drawer.FinishDrawing()
		size = drawer.GetMolSize(scaff_mol)
		##Get the SVG-text
		svg_text_original = drawer.GetDrawingText().strip()
		#Basic SVG-text adjustment
		#New header
		groupID = "<svg xmlns='http://www.w3.org/2000/svg' id = '" + str(scaffold_id_str) + "'  viewBox = '0 0 200 200'  width = '200' height = '200' > <g id = '" + str(scaffold_id_str) + "' opacity = '0.7' class = 'thing'>"
		svg_text = re.sub(r'(?m)<\?xml.*?</rect>', groupID, svg_text_original, flags=re.S)
		if scaff_count[0] == '*': interesting_iresult = svg_text_original
		#Prepare to find the largest molecule out there:
		atom_x_min =  1000000
		atom_x_max = -1000000
		atom_y_min =  1000000
		atom_y_max = -1000000
		#Get the array of strings from svg text
		svg_strings = svg_text.split('\n')
		#Prepare to hold the bonds
		cc_bonds = []
		for bond in scaff_mol.GetBonds():
			bond_id = bond.GetIdx()
			anum1 = bond.GetBeginAtom().GetAtomicNum()
			anum2 = bond.GetEndAtom().GetAtomicNum()
			if anum1 == 6 and anum2 == 6: cc_bonds.append(bond_id)
		typical_bond_length = 0
		#Get xmin, ymin, xmax, ymax for ellipses
		for i, svg_string in enumerate(svg_strings):
			if svg_string.startswith(('<ellipse', '<circle')):
				#Get coordinates and radii
				radius = 0
				radius_x = 0
				radius_y = 0
				elem_x = float(re.search("cx='(.+?)'", svg_string)[1])
				elem_y = float(re.search("cy='(.+?)'", svg_string)[1])
				radius_candidate = re.search("r='(.+?)'", svg_string)
				radius_x_candidate = re.search("rx='(.+?)'", svg_string)
				radius_y_candidate = re.search("ry='(.+?)'", svg_string)
				if radius_candidate:
					radius = float(radius_candidate[1])
					if elem_x + radius > atom_x_max: atom_x_max = elem_x + radius
					if elem_x - radius < atom_x_min: atom_x_min = elem_x - radius
					if elem_y + radius > atom_y_max: atom_y_max = elem_y + radius
					if elem_y - radius < atom_y_min: atom_y_min = elem_y - radius
				if radius_x_candidate and radius_y_candidate:
					radius_x = float(re.search("rx='(.+?)'", svg_string)[1])
					radius_y = float(re.search("ry='(.+?)'", svg_string)[1])
					if elem_x + radius_x > atom_x_max: atom_x_max = elem_x + radius_x
					if elem_x - radius_x < atom_x_min: atom_x_min = elem_x - radius_x
					if elem_y + radius_y > atom_y_max: atom_y_max = elem_y + radius_y
					if elem_y - radius_y < atom_y_min: atom_y_min = elem_y - radius_y
				#Modify elements
				svg_string = svg_string.replace(';stroke:#FFFFFF;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1', ';stroke:none')
			#Get the bond data
			if svg_string.startswith("<path class='bond"):
				#Get the IDs
				bond_id = int(re.search("class='bond-([0-9]+)", svg_string)[1])
				if ( any(cc_bonds == bond_id for cc_bonds in cc_bonds) ):
					#Get the start/end coordinates
					start = re.search("d='M (.+?)'", svg_string)[1]
					if start.count("L")==1 and start.count("Z")==0 :
						end = re.search("L (.+)", start)[1]
						start_x = float( re.search("^[^,]*", start)[0] )
						start_y = float( re.search(",([^,]*) L", start)[1] )
						end_x = float( re.search("L ([^,]*)", start)[1] )
						end_y = float( re.search("L .*,(.*)", start)[1] )
						#Get the distance
						distance = (end_x - start_x)**2  +  (end_y - start_y)**2
						distance_diff_typ1 = abs( 1 - 874/distance )
						distance_diff_typ2 = abs( 1 - 600/distance )
						#Select the most typical bond length
						if typical_bond_length == 0 or min( ( abs(1 - 874/typical_bond_length) ), ( abs(1 - 600/typical_bond_length) ) ) > min( distance_diff_typ1, distance_diff_typ2 ): typical_bond_length = distance
			svg_strings[i] = svg_string
		svg_text = '\n'.join(svg_strings)
		#Calculate the additional score to make the most typical bond length in scaffold equal to the one of the most typical bond lengthes in general (600 and 874)
		if typical_bond_length != 0:
			if ( 1 - 874/typical_bond_length ) < ( 1 - 600/typical_bond_length ):
				add_score = 874 / typical_bond_length
			else:
				add_score = 600 / typical_bond_length
		else:
			add_score = 2
			atom_x_min = 0
			atom_x_max = size[0]
			atom_y_min = 0
			atom_y_max = size[1]
		score = score * add_score
		bbox_scored = [atom_x_min*score, atom_y_min*score, atom_x_max*score, atom_y_max*score]
		scaff_dimension_scored = max(bbox_scored[2] - bbox_scored[0], bbox_scored[3] - bbox_scored[1])
		#Process SVG
		if rings:
			#Get the polygons
			polygon_raw = re.search("<rect.*</svg>", svg_text, re.S).group(0)
			#Delete polygons from the top
			svg_text = svg_text.replace(polygon_raw, "</g>\n</svg>")
			#Clear the polygons
			polygon_raw = re.sub(".*</rect>\n", "", polygon_raw)
			polygon_raw = polygon_raw.replace("</svg>", "")
			#Place the polygons to the bottom
			svg_text = svg_text.replace("'thing'>\n", "'thing'>\n"+polygon_raw)
		else:
			svg_text = svg_text.replace("</svg>", "</g>\n</svg>")
		##Prepare the record to store
		data_to_transform.append({'id': scaffold_id_str, 'scaff_smiles': scaff_count[0], 'structure_ids': scaff_count[2], 'score': score, 'bbox_scored': bbox_scored, 'status': scaffold_status,
							'size': scaff_dimension_scored, 'svg': svg_text})
	##Resize SVGs
	#Sort the SVGs
	data_to_transform = sorted(data_to_transform, key=itemgetter('size'), reverse=True)
	#Transform
	number = 0
	for item in data_to_transform:
		number += 1
		scaffold_id = item['id']
		if number == 1:
			bbox_len = float(item['size'])
			scaling_factor = square_size / bbox_len
		bbox_scaled = [item['bbox_scored'][0]*scaling_factor, item['bbox_scored'][1]*scaling_factor, item['bbox_scored'][2]*scaling_factor, item['bbox_scored'][3]*scaling_factor]
		length_scaled = bbox_scaled[2] - bbox_scaled[0]
		height_scaled = bbox_scaled[3] - bbox_scaled[1]
		transform_factor = (item['score'] * scaling_factor)
		#It is important: "str(-bbox_scaled[0])" to avoid the "--"
		transform_matrix = "(" + str(transform_factor) + ", 0.000000, 0.000000, " + str(transform_factor)+ ", " + str(-bbox_scaled[0]) + ", " + str(-bbox_scaled[1]) + ")"
		svg_text_transform = re.sub(r'(?m)<\?xml.*?</rect>', groupID, item['svg'], flags=re.S)
		#Transform renderable SVG elements
		svg_text_transformed =   svg_text_transform.replace("/>", " transform='matrix" + str(transform_matrix) + "' />")
		svg_text_transformed =   svg_text_transformed.replace("\n", "")
		#Append
		data_to_json.append({'id': item['id'], 'scaff_smiles': item['scaff_smiles'], 'structure_ids': item['structure_ids'], 'min_x': bbox_scaled[0], 'min_y': bbox_scaled[1], 'score': item['score'], 'length': length_scaled, 'height': height_scaled,  'svg_thing': str(base64.b64encode(zlib.compress(svg_text_transformed.encode('utf-8'))))})
	##Gather data for export to browser
	output_json = json.dumps({'meta': meta_data, 'counts': counts, 'depicted_description': selected_description, 'data': data_to_json})
	##Output
	return output_json



try:
	result = func_timeout(1200, process_input, args=(data_input,))
except FunctionTimedOut:
	result = '{"news": "time is out"}'
except Exception as e:
	result = '{"news": "smth is wrong"}'




sys.exit(print(json.dumps(result)))
