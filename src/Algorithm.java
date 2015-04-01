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
		this.pi = 3.1415926535897932384;
		this.final_signatures = new ArrayList<AtomSignature>();
	}
	void data_inizialization(String pathname) { //pathname of a list of sdf file paths with activity to fill data
		//here fill molecule_data, activity_data and allocate memory for atom states and fill it with 0
		//while initializing determine how many activity compounds there are
	}
	class Signature_state {
		Signature_state() {
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
			this.molecule_number = i;
			this.atom_number = j;
		}
	};
	void run() {
		for (int height = 0; height < 15; height++) { //then redo while there are interesting substructures
			Map<AtomSignature, Signature_state> current_signatures= new HashMap<AtomSignature, Signature_state>();
			int molecule_index = 0;
			for (IMolecule m: this.molecule_data) {
				for (int i = 0; i < m.getAtomCount(); i++) {
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
				}
				molecule_index++;
			}
			for (Map.Entry<AtomSignature, Signature_state> substructure: current_signatures.entrySet()) {
				if (substructure.getValue().n_dash > minimum_occurence) {
					double current_frequency = substructure.getValue().m_dash / substructure.getValue().n_dash;
					if (current_frequency > substructure_frequency) {
						double current_p_value = 0;
						double n_ = substructure.getValue().n_dash;
						double m_ = substructure.getValue().m_dash;
						for (double k = m_; k < n_; k++) {
							current_p_value += Math.pow(n_/(2*pi*(n_-k)*k), 0.5)*Math.pow(active_num/k, k)*Math.pow((data_num - active_num)/(n_ - k), n_ - k)*Math.pow(n_/data_num, n_);
						}
						current_p_value += Math.pow(active_num/data_num, n_);
						if (current_p_value < p_value_threshold) {
							for (Signature_position sp: substructure.getValue().met_positions) {
								atom_states.get(sp.molecule_number).set(sp.atom_number, true);
							}
							final_signatures.add(substructure.getKey());
						}
					}
				}
				else {
					for (Signature_position sp: substructure.getValue().met_positions) {
						atom_states.get(sp.molecule_number).set(sp.atom_number, true);
					}
				}
			}
		}
	}
	double p_value_threshold;
	int minimum_occurence;
	double substructure_frequency;
	int data_num;
	int active_num;
	double pi;
	ArrayList<IMolecule> molecule_data;
	ArrayList<Boolean> activity_data;
	ArrayList<ArrayList<Boolean>> atom_states;
	ArrayList<AtomSignature> final_signatures;
}
