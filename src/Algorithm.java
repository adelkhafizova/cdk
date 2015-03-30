import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.signature.AtomSignature;


public class Algorithm {
	Algorithm(double p_value_threshold, int minimum_occurence, double substructure_frequency) {
		this.p_value_threshold = p_value_threshold;
		this.minimum_occurence = minimum_occurence;
		this.substructure_frequency = substructure_frequency;
		this.molecule_data = new ArrayList<IMolecule>();
		this.activity_data = new ArrayList<Boolean>();
		this.atom_states = new ArrayList<ArrayList<Boolean>>(); //true - nothing else to do with it, false - look
	}
	void data_inizialization(String pathname) { //pathname of a list of sdf file paths with activity to fill data
		//here fill molecule_data, activity_data and allocate memory for atom states and fill it with 0
	}
	class Signature_state {
		Signature_state() {
			//this.processed_condition = 0;
			this.met_positions = new ArrayList<Signature_position>();
			n_dash = 0;
			m_dash = 0;
		}
		void add (int i, int j, Boolean activity) {
			Signature_position sp = new Signature_position(i, j);
			this.met_positions.add(sp);
			n_dash++;
			if (activity)
				m_dash++;
		}
		//int processed_condition; //not met; met but not interesting; met but interesting; included in final list
		ArrayList<Signature_position> met_positions; //positions of atoms, where this signature was met
		int n_dash;//number of substructures
		int m_dash;//number of active substructures
	}
	class Signature_position {
		int molecule_number;
		int atom_number;
		Signature_position(int i, int j) {
			this.molecule_number = 0;
			this.atom_number = 0;
		}
	};
	void run() {
		for (int height = 0; height < 15; height++) { //then redo while there are interesting substructures
			Map<AtomSignature, Signature_state> current_signatures= new HashMap<AtomSignature, Signature_state>();
			int molecule_index = 0;
			for (IMolecule m: this.molecule_data) {
				for (int i = 0; i < m.getAtomCount(); i++) {
					//AtomSignature as = new AtomSignature(i, height, m);
					if (this.atom_states.get(molecule_index).get(i) == false) { //not in the current signatures and not viewed
						AtomSignature as = new AtomSignature(i, height, m);
						if (current_signatures.containsKey(as) == false) {
							Signature_state ss = new Signature_state();
							ss.add(molecule_index, i, this.activity_data.get(molecule_index));
							current_signatures.put(as, ss);
						}
						else {
							current_signatures.get(as).add(molecule_index, i, this.activity_data.get(molecule_index));
						}
					}
					/*case true: //met but not interesting
						break;*/
					/*case 2: //met but interesting
						AtomSignature as1 = new AtomSignature(i, height, m);
						if (current_signatures.containsKey(as1) == false) {
							Signature_state ss1 = new Signature_state();
							ss1.add(molecule_index, i, this.activity_data.get(molecule_index));
							current_signatures.put(as1, ss1);
						}
						else {
							current_signatures.get(as1).add(molecule_index, i, this.activity_data.get(molecule_index));
						}
						break;
					case 3: //finished
						break;*/
				}
				molecule_index++;
			}
			
		}
	}
	double p_value_threshold;
	int minimum_occurence;
	double substructure_frequency;
	ArrayList<IMolecule> molecule_data;
	ArrayList<Boolean> activity_data;
	ArrayList<ArrayList<Boolean>> atom_states;
}
