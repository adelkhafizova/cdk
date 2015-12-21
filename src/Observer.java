import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IteratingMDLReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

public class Observer {
    Observer(String original, String trained) {
        original_data_file = original;
        trained_data_file = trained;
    }
    void initialize_algorithms(double p_value_threshold, int minimum_occurrence, double substructure_frequency) {
        algorithm_original = new Algorithm(p_value_threshold, minimum_occurrence, substructure_frequency);
        algorithm_original.data_initialization(original_data_file, "Inhibitor", "Noninhibitor");
        algorithm_trained = new Algorithm(p_value_threshold, minimum_occurrence, substructure_frequency);
        algorithm_trained.data_initialization(trained_data_file, "Inhibitor", "Noninhibitor");
    }
    void analyze() {
        Iterator it = algorithm_original.final_signatures_active.entrySet().iterator();
        while (it.hasNext()) {

        }
        float active = active_active.size()/algorithm_original.final_signatures_active.size();
        float false_negative = active_inactive.size()/algorithm_original.final_signatures_active.size();
        float false_positive = inactive_active.size()/algorithm_original.final_signatures_inactive.size();
        float inactive = inactive_inactive.size()/algorithm_original.final_signatures_inactive.size();
    }
    void run() throws FileNotFoundException, UnsupportedEncodingException, CDKException{
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
