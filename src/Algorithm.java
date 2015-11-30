import java.io.*;
//import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Iterator;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.signature.AtomSignature;
import org.openscience.cdk.signature.MoleculeSignature;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.io.SDFWriter;
//import com.google.common.math.BigIntegerMath;

public class Algorithm {
    Algorithm(double p_value_threshold, int minimum_occurrence, double substructure_frequency) {
        this.p_value_threshold = p_value_threshold;
        this.minimum_occurrence = minimum_occurrence;
        this.substructure_frequency = substructure_frequency;
        this.molecule_data = new ArrayList<IMolecule>();
        this.activity_data = new ArrayList<Boolean>();
        this.atom_states = new ArrayList<ArrayList<Boolean>>(); //true - nothing else to do with it, false - look
        this.pi = 3.1415926535897932384;
        this.final_signatures_active = new HashMap<String, String>();
        this.final_signatures_inactive = new HashMap<String, String>();
        this.heights_active = new HashMap<String, Integer>();
        this.heights_inactive = new HashMap<String, Integer>();
        this.data_num = 0;
        this.active_num = 0;
    }

    void data_initialization(String pathname, String activity_name, String inactivity_name) {
        try {
            File sdfFile = new File(pathname);
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance());
            int inactive_num = 0;
            while (reader.hasNext()) {
                IMolecule molecule = (IMolecule) reader.next();
                molecule_data.add(molecule);
                ArrayList<Boolean> atom_states_local = new ArrayList<Boolean>(molecule.getAtomCount());
                for (int t = 0; t < molecule.getAtomCount(); t++) {
                    atom_states_local.add(false);
                }
                atom_states.add(atom_states_local);
                data_num++;
                String properties = molecule.getProperties().values().toString();
                if (properties.contains(inactivity_name)) {
                    activity_data.add(false);
                    ++inactive_num;
                } else {
                    if (properties.contains(activity_name)) {
                        activity_data.add(true);
                        ++active_num;
                    }
                }
            }
            if (active_num + inactive_num != molecule_data.size()) {
                throw new Exception("Invalid data, not all molecules have activity value");
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

        void add(int i, int j, Boolean activity) {
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
    }

    String atom_signature_to_inchi(AtomSignature as) throws CDKException {
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IAtomContainer ac = MoleculeSignature.fromSignatureString(as.toCanonicalString(), builder);
        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
        InChIGenerator gen = factory.getInChIGenerator(ac);
        String inchi = gen.getInchi();
        return inchi;
    }

    String atom_signature_to_inchi(String as) throws CDKException {
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IAtomContainer ac = MoleculeSignature.fromSignatureString(as, builder);
        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
        InChIGenerator gen = factory.getInChIGenerator(ac);
        String inchi = gen.getInchi();
        return inchi;
    }

    void find_substructures_with_height(Map<String, Signature_state> current_signatures,
                                        Map<String, String> signature_to_inchi, int height) throws CDKException{
        int molecule_index = 0;
        for (IMolecule m : this.molecule_data) { //find all interesting substructures
            for (int i = 0; i < m.getAtomCount(); i++) {
                if (this.atom_states.get(molecule_index).get(i) == false) { //not in the current signatures and not viewed
                    AtomSignature as = new AtomSignature(i, height, m);
                    signature_to_inchi.put(as.toCanonicalString(), atom_signature_to_inchi(as));
                    if (current_signatures.containsKey(as.toString()) == false) {
                        Signature_state ss = new Signature_state();
                        ss.add(molecule_index, i, this.activity_data.get(molecule_index));
                        current_signatures.put(as.toString(), ss);
                    } else {
                        current_signatures.get(as.toString()).add(molecule_index, i, this.activity_data.get(molecule_index));
                    }
                }
            }
            molecule_index++;
        }
    }

    void process_signatures(Map<String, Signature_state> current_signatures, int height) throws CDKException {
        for (Map.Entry<String, Signature_state> substructure : current_signatures.entrySet()) {
            if (substructure.getValue().n_dash > minimum_occurrence) {
                double current_frequency = substructure.getValue().m_dash / substructure.getValue().n_dash;
                double n_ = substructure.getValue().n_dash;
                double m_ = substructure.getValue().m_dash;
                if (current_frequency > substructure_frequency) {
                    double current_p_value_active = 0;
                    for (double k = m_; k < n_; k++) {
                        current_p_value_active += Math.pow(n_ / (2 * pi * (n_ - k) * k), 0.5) * Math.pow(active_num / k, k) *
                                Math.pow((data_num - active_num) / (n_ - k), n_ - k) * Math.pow(n_ / data_num, n_);
                        //BigInteger coef = BigIntegerMath.binomial((int) n_, (int) k);
                        //true_p_value += (coef.doubleValue() * Math.pow(active_num, k)) / Math.pow(data_num, k) *
                        //        Math.pow((data_num - active_num), n_ - k) / Math.pow(data_num, n_ - k);
                    }
                    current_p_value_active += Math.pow(active_num / data_num, n_);
                    if (current_p_value_active < p_value_threshold) {
                        for (Signature_position sp : substructure.getValue().met_positions) {
                            atom_states.get(sp.molecule_number).set(sp.atom_number, true);
                        }

                        String as_string = substructure.getKey();
                        String InChi = atom_signature_to_inchi(as_string);
                        final_signatures_active.put(as_string, InChi);
                        heights_active.put(as_string, height);
                    }
                }

                current_frequency = (substructure.getValue().n_dash - substructure.getValue().m_dash) / substructure.getValue().n_dash;

                if (current_frequency > substructure_frequency) {
                    double current_p_value_inactive = 0;
                    m_ = n_ - substructure.getValue().m_dash;
                    for (double k = m_; k < n_; k++) {
                        current_p_value_inactive += Math.pow(n_ / (2 * pi * (n_ - k) * k), 0.5) * Math.pow((data_num - active_num) / k, k) *
                                Math.pow(active_num / (n_ - k), n_ - k) * Math.pow(n_ / data_num, n_);
                    }
                    current_p_value_inactive += Math.pow((data_num - active_num) / data_num, n_);
                    if (current_p_value_inactive < p_value_threshold) {
                        for (Signature_position sp : substructure.getValue().met_positions) {
                            atom_states.get(sp.molecule_number).set(sp.atom_number, true);
                        }
                        String as_string = substructure.getKey();
                        String InChi = atom_signature_to_inchi(as_string);
                        final_signatures_inactive.put(as_string, InChi);
                        heights_inactive.put(as_string, height);
                    }
                }
            } else {
                for (Signature_position sp : substructure.getValue().met_positions) {
                    atom_states.get(sp.molecule_number).set(sp.atom_number, true);
                }
            }
        }
    }

    void dump_data(String path) throws FileNotFoundException, UnsupportedEncodingException {
        PrintWriter writer = new PrintWriter(path + "final_signatures_active.txt", "UTF-8");
        for (Map.Entry<String, String> substructure : final_signatures_active.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
        writer = new PrintWriter(path + "final_signatures_inactive.txt", "UTF-8");
        for (Map.Entry<String, String> substructure : final_signatures_inactive.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
        writer = new PrintWriter(path + "signatures_active_heights.txt", "UTF-8");
        for (Map.Entry<String, Integer> substructure : heights_active.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
        writer = new PrintWriter(path + "signatures_inactive_heights.txt", "UTF-8");
        for (Map.Entry<String, Integer> substructure : heights_inactive.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
    }

    void run(String path) throws FileNotFoundException, UnsupportedEncodingException, CDKException {
        for (int height = 0; height < 15; height++) { //then redo while there are interesting substructures
            Map<String, Signature_state> current_signatures = new HashMap<String, Signature_state>();
            Map<String, String> signature_to_inchi = new HashMap<String, String>();
            find_substructures_with_height(current_signatures, signature_to_inchi, height);
            process_signatures(current_signatures, height);
        }
        dump_data(path);
    }

	double p_value_threshold;
	int minimum_occurrence;
	double substructure_frequency;
	int data_num;
	int active_num;
	double pi;
	ArrayList<IMolecule> molecule_data;
	ArrayList<Boolean> activity_data;
	ArrayList<ArrayList<Boolean>> atom_states;
	HashMap<String, String> final_signatures_active;
	HashMap<String, String> final_signatures_inactive;
    HashMap<String, Integer> heights_active;
    HashMap<String, Integer> heights_inactive;
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException, CDKException,
            IOException {
	Algorithm a = new Algorithm(0.05, 5, 0.8);
	a.data_initialization("cas_4337.sdf", "mutagen", "nonmutagen");
	a.run(path);
        Classification classifier = new Classification();
        classifier.data_initialization("");
        classifier.classify_data("test_mol_for_class.sdf");
        //a.classify_data("final_signatures_active.txt", "final_signatures_inactive.txt", "test_mol_for_class.sdf");
	}
}
