import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.SystemOutLoggingTool;

import java.io.*;
import java.util.*;

public class Observer {
    Observer(String original, String trained) {
        original_data_file = original;
        trained_data_file = trained;
    }
    void initialize_algorithms(double p_value_threshold, int minimum_occurrence, double substructure_frequency) {
        algorithm_original = new Algorithm(p_value_threshold, minimum_occurrence, substructure_frequency,
                                           new HashMap<String, Double>(), new HashMap<String, Double>());
        algorithm_original.data_initialization(original_data_file, "active", "inactive");
        algorithm_trained = new Algorithm(p_value_threshold, minimum_occurrence, substructure_frequency,
                                          algorithm_original.final_signatures_active,
                                          algorithm_original.final_signatures_inactive);
        algorithm_trained.data_initialization(trained_data_file, "active", "inactive");
    }
    void analyze() throws FileNotFoundException, UnsupportedEncodingException, CDKException {
        active_active = new HashMap<String, String>();
        active_inactive = new HashMap<String, String>();
        inactive_active = new HashMap<String, String>();
        inactive_inactive = new HashMap<String, String>();
        HashSet<String> active_remained_active = new HashSet<String>();
        HashSet<String> active_became_inactive = new HashSet<String>();
        HashSet<String> active_became_insignificant = new HashSet<String>();
        HashSet<String> insignificant_became_active = new HashSet<String>();
        HashSet<String> inactive_remained_inactive = new HashSet<String>();
        HashSet<String> inactive_became_active = new HashSet<String>();
        HashSet<String> insignificant_became_inactive = new HashSet<String>();
        Iterator it = algorithm_trained.active.entrySet().iterator();
        HashMap<String, Integer> plus_plus = new HashMap<String, Integer>();
        HashMap<String, Integer> plus_minus = new HashMap<String, Integer>();
        HashMap<String, Integer> minus_minus = new HashMap<String, Integer>();
        HashMap<String, Integer> minus_minus = new HashMap<String, Integer>();
        float active = 0;
        float false_negative = 0;
        float false_positive = 0;
        float inactive = 0;
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            if (algorithm_original.active.containsKey(pair.getKey())) {
                active_active.put((String) pair.getKey(), (String) pair.getValue());
                ++active;
            }
            if (algorithm_original.inactive.containsKey(pair.getKey())) {
                active_inactive.put((String)pair.getKey(), (String)pair.getValue());
                ++false_positive;
            }
        }
        it = algorithm_trained.inactive.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            if (algorithm_original.active.containsKey(pair.getKey())) {
                inactive_active.put((String)pair.getKey(), (String)pair.getValue());
                ++false_negative;
            }
            if (algorithm_original.inactive.containsKey(pair.getKey())) {
                inactive_inactive.put((String)pair.getKey(), (String)pair.getValue());
                ++inactive;
            }
        }
        it = algorithm_original.active.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            if (algorithm_trained.active.containsKey(pair.getKey())) {
                active_remained_active.add((String) pair.getKey());
            }
            else {
                if (algorithm_trained.inactive.containsKey(pair.getKey())) {
                    active_became_inactive.add((String) pair.getKey());
                } else {
                    active_became_insignificant.add((String)pair.getKey());
                }
            }
        }
        it = algorithm_original.inactive.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            if (algorithm_trained.inactive.containsKey(pair.getKey())) {
                inactive_remained_inactive.add((String) pair.getKey());
            }
            else {
                if (algorithm_trained.active.containsKey(pair.getKey())) {
                    inactive_became_active.add((String) pair.getKey());
                } else {
                    inactive_became_insignificant.add((String)pair.getKey());
                }
            }
        }
        it = algorithm_trained.active.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            if (algorithm_original.active.containsKey(pair.getKey()) == False and algorithm_original.inactive.containsKey(pair.getKey())) {
                insignificant_became_active.add((String) pair.getKey());
            }
        }
        it = algorithm_trained.inactive.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            if (algorithm_original.active.containsKey(pair.getKey()) == False and algorithm_original.inactive.containsKey(pair.getKey())) {
                insignificant_became_inactive.add((String) pair.getKey());
            }
        }
        PrintWriter writer = new PrintWriter("categorized_signatures.txt", "UTF-8");
        print_set(active_remained_active, 0, writer);
        print_set(active_became_inactive, 1, writer);
        print_set(active_became_insignificant, 2, writer);
        print_set(inactive_remained_inactive, 3, writer);
        print_set(inactive_became_active, 4, writer);
        print_set(inactive_became_insignificant, 5, writer);
        print_set(insignificant_became_active, 6, writer);
        print_set(insignificant_became_inactive, 7, writer);
        writer.close();
        //float active = (float)active_active.size()/(float)algorithm_original.final_signatures_active.size();
        //float false_negative = (float)active_inactive.size()/(float)algorithm_original.final_signatures_active.size();
        //float false_positive = (float)inactive_active.size()/(float)algorithm_original.final_signatures_inactive.size();
        //float inactive = (float)inactive_inactive.size()/(float)algorithm_original.final_signatures_inactive.size();
        System.out.println(active);
        System.out.println(false_negative);
        System.out.println(false_positive);
        System.out.println(inactive);
    }
    void print_set(HashSet<String> molecules, Integer category, PrintWriter writer) throws CDKException {
        Iterator<String> it = molecules.iterator();
        while (it.hasNext()) {
            InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
            InChIToStructure intostruct = factory.getInChIToStructure(
                    it.next().toString(), DefaultChemObjectBuilder.getInstance()
            );

            IAtomContainer container = intostruct.getAtomContainer();

            IMolecule molecule = new Molecule(container);

            SmilesGenerator sg = new SmilesGenerator();
            String smiles = sg.createSMILES(molecule);
            writer.println(smiles + '\t' + category.toString());
        }
    }
    void run() throws Exception {
        algorithm_original.run("");
        algorithm_trained.run("");
    }
    Algorithm algorithm_original;
    Algorithm algorithm_trained;
    String original_data_file;
    String trained_data_file;
    HashMap<String, String> active_active; //active originally and active in trained
    HashMap<String, String> active_inactive; //active originally but inactive in trained
    HashMap<String, String> inactive_inactive; // inactive originally and inactive in trained
    HashMap<String, String> inactive_active; //inactive originally and active in trained
}
