import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Iterator;

import com.sun.org.apache.xpath.internal.operations.Bool;
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

import com.google.common.math.BigIntegerMath;

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

    void data_initialization(String pathname, String inactivity_name) {
        try {
            File sdfFile = new File(pathname);
            IteratingMDLReader reader = new IteratingMDLReader(new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance());
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
                } else {
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

    void run() throws FileNotFoundException, UnsupportedEncodingException, CDKException {
        for (int height = 0; height < 15; height++) { //then redo while there are interesting substructures
            Map<String, Signature_state> current_signatures = new HashMap<String, Signature_state>();
            Map<String, String> signature_to_inchi = new HashMap<String, String>();
            int molecule_index = 0;
            for (IMolecule m : this.molecule_data) {
                for (int i = 0; i < m.getAtomCount(); i++) {
                    if (this.atom_states.get(molecule_index).get(i) == false) { //not in the current signatures and not viewed
                        AtomSignature as = new AtomSignature(i, height, m);
                        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
                        IAtomContainer ac = MoleculeSignature.fromSignatureString(as.toCanonicalString(), builder);
                        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
                        InChIGenerator gen = factory.getInChIGenerator(ac);
                        String inchi = gen.getInchi();
                        signature_to_inchi.put(as.toCanonicalString(), inchi);
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

            for (Map.Entry<String, Signature_state> substructure : current_signatures.entrySet()) {
                if (substructure.getValue().n_dash > minimum_occurrence) {
                    double current_frequency = substructure.getValue().m_dash / substructure.getValue().n_dash;
                    double n_ = substructure.getValue().n_dash;
                    double m_ = substructure.getValue().m_dash;
                    if (current_frequency > substructure_frequency) {
                        double current_p_value_active = 0;
                        //double true_p_value = 0;
                        for (double k = m_; k < n_; k++) {
                            current_p_value_active += Math.pow(n_ / (2 * pi * (n_ - k) * k), 0.5) * Math.pow(active_num / k, k) *
                                    Math.pow((data_num - active_num) / (n_ - k), n_ - k) * Math.pow(n_ / data_num, n_);
                            //BigInteger coef = BigIntegerMath.binomial((int) n_, (int) k);
                            //true_p_value += (coef.doubleValue() * Math.pow(active_num, k)) / Math.pow(data_num, k) *
                            //        Math.pow((data_num - active_num), n_ - k) / Math.pow(data_num, n_ - k);
                        }
                        current_p_value_active += Math.pow(active_num / data_num, n_);
                        //true_p_value += Math.pow(active_num / data_num, n_);
                        //System.out.println(current_p_value_active);
                        //System.out.println(true_p_value);
                        //System.out.println("next");
                        if (current_p_value_active < p_value_threshold) {
                            for (Signature_position sp : substructure.getValue().met_positions) {
                                atom_states.get(sp.molecule_number).set(sp.atom_number, true);
                            }
                            IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
                            String as_string = substructure.getKey();
                            IAtomContainer ac = MoleculeSignature.fromSignatureString(as_string, builder);
                            InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
                            InChIGenerator gen = factory.getInChIGenerator(ac);
                            String InChi = gen.getInchi();
                            final_signatures_active.put(as_string, InChi);
                            heights_active.put(as_string, height);
                            /*if (substructure.getKey().equals("[H]([C]([C]([H][H][N])[Cl][H]))")) {
                                //System.out.println(true_p_value);
                                System.out.println(n_);
                                System.out.println(m_);
                            }*/
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
                            IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
                            String as_string = substructure.getKey();
                            IAtomContainer ac = MoleculeSignature.fromSignatureString(as_string, builder);
                            InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
                            InChIGenerator gen = factory.getInChIGenerator(ac);
                            String InChi = gen.getInchi();
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
        System.out.println(final_signatures_active.size());
        System.out.println(final_signatures_active.get("[Cl]([C]([C]([H][H][N])[H][H]))"));
        PrintWriter writer = new PrintWriter("final_signatures_active.txt", "UTF-8");
        for (Map.Entry<String, String> substructure : final_signatures_active.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
        writer = new PrintWriter("final_signatures_inactive.txt", "UTF-8");
        for (Map.Entry<String, String> substructure : final_signatures_inactive.entrySet()) {
            writer.println(substructure.getKey() + ' ' + substructure.getValue());
        }
        writer.close();
    }

    void classify_data() throws CDKException, java.io.FileNotFoundException, java.io.UnsupportedEncodingException  {
        ArrayList<Boolean> active = new ArrayList<Boolean>();
        ArrayList<Boolean> inactive = new ArrayList<Boolean>();
        ArrayList<HashSet<String>> active_fragments = new ArrayList<HashSet<String>>();
        ArrayList<HashSet<String> > inactive_fragments = new ArrayList<HashSet<String>>();
        int index = 0;
        System.out.println("Progress bar");
        for (IMolecule m : this.molecule_data) {
            active.add(false);
            inactive.add(false);
            active_fragments.add(new HashSet<String>());
            inactive_fragments.add(new HashSet<String>());
            for (int i = 0; i < m.getAtomCount(); ++i) {
                for (Map.Entry<String, String> substructure: final_signatures_active.entrySet()) {
                    String as = substructure.getKey();
                    int begin = as.indexOf('[');
                    int end = as.indexOf(']', begin + 1);
                    if (m.getAtom(i).getSymbol().equals(as.substring(begin + 1, end))) {
                        AtomSignature atom_s = new AtomSignature(i, heights_active.get(as), m);
                        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
                        String as_string = atom_s.toString();
                        IAtomContainer ac = MoleculeSignature.fromSignatureString(as_string, builder);
                        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
                        InChIGenerator gen = factory.getInChIGenerator(ac);
                        String InChi = gen.getInchi();
                        if (InChi.equals(substructure.getValue())) {
                            active_fragments.get(index).add(InChi);
                            //System.out.println("Found active");
                            active.set(index, true);
                            break;
                        }
                    }
                }
                for (Map.Entry<String, String> substructure: final_signatures_inactive.entrySet()) {
                    String as = substructure.getKey();
                    int begin = as.indexOf('[');
                    int end = as.indexOf(']', begin + 1);
                    if (m.getAtom(i).getSymbol().equals(as.substring(begin + 1, end))) {
                        AtomSignature atom_s = new AtomSignature(i, heights_inactive.get(as), m);
                        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
                        String as_string = atom_s.toString();
                        IAtomContainer ac = MoleculeSignature.fromSignatureString(as_string, builder);
                        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
                        InChIGenerator gen = factory.getInChIGenerator(ac);
                        String InChi = gen.getInchi();
                        if (InChi.equals(substructure.getValue())) {
                            inactive_fragments.get(index).add(InChi);
                            //System.out.println("Found inactive");
                            inactive.set(index, true);
                            break;
                        }
                    }
                }
            }
            ++index;
            System.out.println(((double) index) / this.molecule_data.size() * 100);
        }
        int active_count = 0;
        int inactive_count = 0;
        ArrayList<Integer> final_activity = new ArrayList<Integer>();
        for (int i = 0; i < active.size(); ++i) {
            if (active.get(i)) {
                if (inactive.get(i)) {
                    System.out.println("unrecognized");
                    final_activity.add(2);
                } else {
                    System.out.println("mutagen");
                    ++active_count;
                    final_activity.add(1);
                }
            } else {
                System.out.println("nonmutagen");
                ++inactive_count;
                final_activity.add(0);
            }
        }
        System.out.println("Inactive");
        for (Boolean activity: inactive) {
            System.out.println(activity);
        }
        System.out.println("Active count" + active_count);
        System.out.println("Inactive count" + inactive_count);
        PrintWriter output_file = new PrintWriter("final_results", "UTF-8");
        for (int i = 0; i < molecule_data.size(); ++i) {
            output_file.println(molecule_data.get(i).toString() + '\n');
            if (final_activity.get(i) == 0) {
                output_file.println("nonmutagen" + '\n');
            } else {
                if (final_activity.get(i) == 1)  {
                    output_file.println("mutagen" + '\n');
                } else {
                    output_file.println("unrecognized" + '\n');
                }
            }
            Iterator<String> itr = active_fragments.get(i).iterator();
            while (itr.hasNext()) {
                output_file.print(itr.next() + ' ');
            }
            output_file.println();
            itr = inactive_fragments.get(i).iterator();
            while (itr.hasNext()) {
                output_file.print(itr.next() + ' ');
            }
            output_file.println();
        }
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
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException, CDKException {
		Algorithm a = new Algorithm(0.05, 5, 0.8);
		a.data_initialization("cas_4337.sdf", "nonmutagen");
		a.run();
        a.classify_data();
	}
}
