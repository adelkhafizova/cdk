import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class Observer {
    Observer(String original, String trained, String activity_name, String inactivity_name) {
        original_data_file = original;
        trained_data_file = trained;
        this.activity_name = activity_name;
        this.inactivity_name = inactivity_name;
    }
    void initialize_algorithms(double p_value_threshold, int minimum_occurrence, double substructure_frequency, boolean use_accuracies) throws Exception {
        algorithm_original = new Algorithm(p_value_threshold, minimum_occurrence, substructure_frequency,
                new HashMap<String, Double>(), new HashMap<String, Double>(), true);
        algorithm_original.data_initialization(original_data_file, activity_name, inactivity_name, use_accuracies);
        algorithm_trained = new Algorithm(p_value_threshold, minimum_occurrence, substructure_frequency,
                algorithm_original.final_signatures_active,
                algorithm_original.final_signatures_inactive, false);
        algorithm_trained.data_initialization(trained_data_file, activity_name, inactivity_name, use_accuracies);
    }

    void analyze(Map<String, Algorithm.SignatureInfo> table) throws IOException, CDKException {
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
        HashSet<String> inactive_became_insignificant = new HashSet<String>();
        HashSet<String> insignificant_became_inactive = new HashSet<String>();
        Iterator it = algorithm_trained.active.entrySet().iterator();
        HashMap<String, Integer> plus_plus = new HashMap<String, Integer>();
        HashMap<String, Integer> plus_minus = new HashMap<String, Integer>();
        HashMap<String, Integer> minus_minus = new HashMap<String, Integer>();
        HashMap<String, Integer> minus_plus = new HashMap<String, Integer>();
        float active = 0;
        float false_negative = 0;
        float false_positive = 0;
        float inactive = 0;
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry) it.next();
            if (algorithm_original.active.containsKey(pair.getKey())) {
                active_active.put((String) pair.getKey(), (String) pair.getValue());
                ++active;
            }
            if (algorithm_original.inactive.containsKey(pair.getKey())) {
                active_inactive.put((String) pair.getKey(), (String) pair.getValue());
                ++false_positive;
            }
        }
        it = algorithm_trained.inactive.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry) it.next();
            if (algorithm_original.active.containsKey(pair.getKey())) {
                inactive_active.put((String) pair.getKey(), (String) pair.getValue());
                ++false_negative;
            }
            if (algorithm_original.inactive.containsKey(pair.getKey())) {
                inactive_inactive.put((String) pair.getKey(), (String) pair.getValue());
                ++inactive;
            }
        }
        it = algorithm_original.active.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry) it.next();
            if (algorithm_trained.active.containsKey(pair.getKey())) {
                active_remained_active.add((String) pair.getKey());
            } else {
                if (algorithm_trained.inactive.containsKey(pair.getKey())) {
                    active_became_inactive.add((String) pair.getKey());
                } else {
                    active_became_insignificant.add((String) pair.getKey());
                }
            }
        }
        it = algorithm_original.inactive.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry) it.next();
            if (algorithm_trained.inactive.containsKey(pair.getKey())) {
                inactive_remained_inactive.add((String) pair.getKey());
            } else {
                if (algorithm_trained.active.containsKey(pair.getKey())) {
                    inactive_became_active.add((String) pair.getKey());
                } else {
                    inactive_became_insignificant.add((String) pair.getKey());
                }
            }
        }
        it = algorithm_trained.active.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry) it.next();
            if (algorithm_original.active.containsKey(pair.getKey()) == false && algorithm_original.inactive.containsKey(pair.getKey()) == false){
                insignificant_became_active.add((String) pair.getKey());
            }
        }
        it = algorithm_trained.inactive.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry) it.next();
            if (algorithm_original.active.containsKey(pair.getKey()) == false && algorithm_original.inactive.containsKey(pair.getKey()) == false){
                insignificant_became_inactive.add((String) pair.getKey());
            }
        }
        PrintWriter writer = new PrintWriter("categorized_signatures.txt", "UTF-8");
        File output = new File("categories.sdf");
        SDFWriter sdf_writer = new SDFWriter(new FileWriter(output));
        print_set(active_remained_active, 0, writer, sdf_writer, table);
        print_set(active_became_inactive, 1, writer, sdf_writer, table);
        print_set(active_became_insignificant, 2, writer, sdf_writer, table);
        print_set(inactive_remained_inactive, 3, writer, sdf_writer, table);
        print_set(inactive_became_active, 4, writer, sdf_writer, table);
        print_set(inactive_became_insignificant, 5, writer, sdf_writer, table);
        print_set(insignificant_became_active, 6, writer, sdf_writer, table);
        print_set(insignificant_became_inactive, 7, writer, sdf_writer, table);
        print_table(table);
        writer.close();
        sdf_writer.close();
        //float active = (float)active_active.size()/(float)algorithm_original.final_signatures_active.size();
        //float false_negative = (float)active_inactive.size()/(float)algorithm_original.final_signatures_active.size();
        //float false_positive = (float)inactive_active.size()/(float)algorithm_original.final_signatures_inactive.size();
        //float inactive = (float)inactive_inactive.size()/(float)algorithm_original.final_signatures_inactive.size();
        System.out.println("Active active size" + active);
        System.out.println("Inactive became active" + false_negative);
        System.out.println("Active became inactive" + false_positive);
        System.out.println("Inactive became inactive" + inactive);
    }

    void print_set(HashSet<String> molecules, Integer category, PrintWriter writer, SDFWriter sdfwriter, Map<String, Algorithm.SignatureInfo> table) throws CDKException {
        Iterator<String> it = molecules.iterator();
        Integer counter = 0;
        if (molecules.size() == 0) {
            return;
        }
        System.out.println(molecules.size());
        while (it.hasNext()) {
            InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
            String next_string = it.next();
            if (next_string == null) {
                System.err.println("Null here!");
                System.err.println("Counter " + counter.toString());
                continue;
            }
            String current_inchi = next_string.toString();
            table.get(current_inchi).fill_category(category);
            InChIToStructure intostruct = factory.getInChIToStructure(
                    current_inchi, DefaultChemObjectBuilder.getInstance()
            );

            IAtomContainer container = intostruct.getAtomContainer();

            IMolecule molecule = new Molecule(container);
            molecule.setProperty("Category", category.toString());
            sdfwriter.write(molecule);
            SmilesGenerator sg = new SmilesGenerator();
            String smiles = sg.createSMILES(molecule);
            if (table.get(current_inchi) == null) {
                System.err.println("No current inchi in table");
            }
            table.get(current_inchi).fill_smiles(smiles);
            writer.println(smiles + '\t' + category.toString());
        }
    }

    void print_table(Map<String, Algorithm.SignatureInfo> table) throws FileNotFoundException, UnsupportedEncodingException {
        Iterator it = table.entrySet().iterator();
        PrintWriter writer = new PrintWriter("statistics.csv", "UTF-8");
        writer.print("Smiles");
        writer.print("\t");
        writer.print("p_value_original");
        writer.print("\t");
        writer.print("p_value_trained");
        writer.print("\t");
        writer.print("active_entry_original");
        writer.print("\t");
        writer.print("inactive_entry_original");
        writer.print("\t");
        writer.print("active_entry_trained");
        writer.print("\t");
        writer.print("inactive_entry_trained");
        writer.print("\t");
        writer.print("category");
        writer.print("\t");
        writer.print("accuracy");
        writer.print("\n");
        DecimalFormat numberFormat = new DecimalFormat("#.00");
        //System.out.println(numberFormat.format(number));
        while (it.hasNext()) {
            Map.Entry<String, Algorithm.SignatureInfo> next_entry = (Map.Entry<String, Algorithm.SignatureInfo>)it.next();
            writer.print(next_entry.getValue().smiles);
            writer.print("\t");
            writer.print(formatNum(next_entry.getValue().p_value_original));
            writer.print("\t");
            writer.print(formatNum(next_entry.getValue().p_value_trained));
            writer.print("\t");
            writer.print(next_entry.getValue().active_entry_count_original);
            writer.print("\t");
            writer.print(next_entry.getValue().inactive_entry_count_original);
            writer.print("\t");
            writer.print(next_entry.getValue().active_entry_count_trained);
            writer.print("\t");
            writer.print(next_entry.getValue().inactive_entry_count_trained);
            writer.print("\t");
            writer.print(next_entry.getValue().category);
            writer.print("\t");
            writer.print(next_entry.getValue().average_accuracy);
            writer.print("\n");
        }

    }

    private static final int MAX_LENGTH = 7;

    private static String formatNum(double number) {
    int digitsAvailable = MAX_LENGTH - 2;
    if (Math.abs(number) < Math.pow(10, digitsAvailable)
            && Math.abs(number) > Math.pow(10, -digitsAvailable)) {
        String format = "0.";
        double temp = number;
        for (int i = 0; i < digitsAvailable; i++) {
            if ((temp /= 10) < 1) {
                format += "#";
            }
        }
        return new DecimalFormat(format).format(number);
    }
    String format = "0.";
    for (int i = 0; i < digitsAvailable; i++) {
            format += "#";
    }
    String r = new DecimalFormat(format + "E0").format(number);
    int lastLength = r.length() + 1;
    while (r.length() > MAX_LENGTH && lastLength > r.length()) {
        lastLength = r.length();
        r = r.replaceAll("\\.?[0-9]E", "E");
    }
    return r;
}



    void run(Map<String, Algorithm.SignatureInfo> table) throws Exception {
        algorithm_original.run(table, "");
        algorithm_trained.run(table, "");
    }
    Algorithm algorithm_original;
    Algorithm algorithm_trained;
    String original_data_file;
    String trained_data_file;
    String activity_name;
    String inactivity_name;
    HashMap<String, String> active_active; //active originally and active in trained
    HashMap<String, String> active_inactive; //active originally but inactive in trained
    HashMap<String, String> inactive_inactive; // inactive originally and inactive in trained
    HashMap<String, String> inactive_active; //inactive originally and active in trained
}
