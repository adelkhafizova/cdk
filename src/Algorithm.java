import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.signature.AtomSignature;

import com.google.common.math.BigIntegerMath;
import com.google.common.math.DoubleMath;

public class Algorithm {
	Algorithm(double p_value_threshold, int minimum_occurence, double substructure_frequency) {
		this.p_value_threshold = p_value_threshold;
		this.minimum_occurence = minimum_occurence;
		this.substructure_frequency = substructure_frequency;
		this.molecule_data = new ArrayList<IMolecule>();
		this.activity_data = new ArrayList<Boolean>();
		this.atom_states = new ArrayList<ArrayList<Boolean>>(); //true - nothing else to do with it, false - look
		this.pi = 3.1415926535897932384;
		this.final_signatures_active = new ArrayList<String>();
		this.final_signatures_inactive = new ArrayList<String>();
		this.data_num = 0;
		this.active_num = 0;
	}
	void data_inizialization(String pathname, String inactivity_name) { 
        try {
        	 File sdfFile = new File(pathname);
        	 IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance());
        	 while (reader.hasNext()) {
        		 IMolecule molecule = (IMolecule)reader.next();
        		 molecule_data.add(molecule);
        		 ArrayList<Boolean> atom_states_local =  new ArrayList<Boolean>(molecule.getAtomCount());
        		 for (int t = 0; t < molecule.getAtomCount(); t++) {
        			 atom_states_local.add(false);
        		 }
        		 atom_states.add(atom_states_local);
        		 data_num++;
        		 String properties = molecule.getProperties().values().toString();
        		 if (properties.indexOf(inactivity_name) != -1) {
        			 activity_data.add(false);
        		 }
        		 else {
        			 activity_data.add(true);
        			 active_num++;
        		 }
        	}
		//here fill molecule_data, activity_data and allocate memory for atom states and fill it with 0
		//while initializing determine how many activity compounds there are
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
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
		double n_dash;//number of substructures
		double m_dash;//number of active substructures
	}
	class Signature_position {
		int molecule_number;
		int atom_number;
		Signature_position(int i, int j) {
			this.molecule_number = i;
			this.atom_number = j;
		}
	};
	void run() throws FileNotFoundException, UnsupportedEncodingException {
		for (int height = 0; height < 15; height++) { //then redo while there are interesting substructures
			Map<String, Signature_state> current_signatures= new HashMap<String, Signature_state>();
			int molecule_index = 0;
			for (IMolecule m: this.molecule_data) {
				for (int i = 0; i < m.getAtomCount(); i++) {
					if (this.atom_states.get(molecule_index).get(i) == false) { //not in the current signatures and not viewed
						AtomSignature as = new AtomSignature(i, height, m);
						if (current_signatures.containsKey(as.toString()) == false) {
							Signature_state ss = new Signature_state();
							ss.add(molecule_index, i, this.activity_data.get(molecule_index));
							current_signatures.put(as.toString(), ss);
						}
						else {
							current_signatures.get(as.toString()).add(molecule_index, i, this.activity_data.get(molecule_index));
						}
					}
				}
				molecule_index++;
			}
				//System.out.println(current_signatures.containsKey("[N]([C]([C]([Cl])))"));
			//System.out.println(current_signatures.size());
			for (Map.Entry<String, Signature_state> substructure: current_signatures.entrySet()) {
				//System.out.println(substructure.getValue().n_dash);
				if (substructure.getValue().n_dash > minimum_occurence) {
					double current_frequency = substructure.getValue().m_dash / substructure.getValue().n_dash;
					double n_ = substructure.getValue().n_dash;
					double m_ = substructure.getValue().m_dash;
					if (current_frequency > substructure_frequency) {
						double current_p_value_active = 0;
						double true_p_value = 0;
						for (double k = m_; k < n_; k++) {
							current_p_value_active += Math.pow(n_/(2*pi*(n_-k)*k), 0.5)*Math.pow(active_num/k, k)*
									Math.pow((data_num - active_num)/(n_ - k), n_ - k)*Math.pow(n_/data_num, n_);
							BigInteger coef = BigIntegerMath.binomial((int) n_, (int) k);
							true_p_value += (coef.doubleValue()*Math.pow(active_num, k))/Math.pow(data_num, k)*
									Math.pow((data_num - active_num), n_ - k)/Math.pow(data_num, n_ - k);
//							System.out.println((coef.doubleValue()*Math.pow(active_num, k))/Math.pow(data_num, k)*
//									Math.pow((data_num - active_num), n_ - k)/Math.pow(data_num, n_ - k));
						}
						current_p_value_active += Math.pow(active_num/data_num, n_);
						true_p_value += Math.pow(active_num/data_num, n_);
						//System.out.println(current_p_value_active);
						//System.out.println(true_p_value);
						//System.out.println("next");
						if (current_p_value_active < p_value_threshold) {
							for (Signature_position sp: substructure.getValue().met_positions) {
								atom_states.get(sp.molecule_number).set(sp.atom_number, true);
							}
							final_signatures_active.add(substructure.getKey());
							if (substructure.getKey().equals("[H]([C]([C]([H][H][N])[Cl][H]))")) {
								System.out.println(true_p_value);
								System.out.println(n_);
								System.out.println(m_);
							}
						}
					}
					
					current_frequency = (substructure.getValue().n_dash - substructure.getValue().m_dash) / substructure.getValue().n_dash;
						
					if (current_frequency > substructure_frequency) {
						double current_p_value_inactive = 0;
						m_ = n_ - substructure.getValue().m_dash;
						for (double k = m_; k < n_; k++) {
								current_p_value_inactive += Math.pow(n_/(2*pi*(n_-k)*k), 0.5)*Math.pow((data_num - active_num)/k, k)*
									Math.pow(active_num/(n_ - k), n_ - k)*Math.pow(n_/data_num, n_);
						}
						current_p_value_inactive += Math.pow((data_num - active_num)/data_num, n_);
						//System.out.println(current_p_value);
						if (current_p_value_inactive < p_value_threshold) {
							for (Signature_position sp: substructure.getValue().met_positions) {
								atom_states.get(sp.molecule_number).set(sp.atom_number, true);
							}
							final_signatures_inactive.add(substructure.getKey());
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
		System.out.println(final_signatures_active.size());
		System.out.println(final_signatures_active.contains("[Cl]([C]([C]([H][H][N])[H][H]))"));
		PrintWriter writer = new PrintWriter("final_signatures_active.txt", "UTF-8");
		for (int w = 0; w < final_signatures_active.size(); w++) {
			writer.println(final_signatures_active.get(w));
		}
		writer.close();
		writer = new PrintWriter("final_signatures_inactive.txt", "UTF-8");
		for (int w = 0; w < final_signatures_inactive.size(); w++) {
			writer.println(final_signatures_inactive.get(w));
		}
		writer.close();
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
	ArrayList<String> final_signatures_active;
	ArrayList<String> final_signatures_inactive;
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		Algorithm a = new Algorithm(0.1, 5, 0.8);
		a.data_inizialization("cas_4337.sdf", "nonmutagen");
		a.run();
	}
}
